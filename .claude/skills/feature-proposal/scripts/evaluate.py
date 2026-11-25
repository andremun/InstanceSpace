"""Run the feature-proposal checks over a set of candidate features.

    python evaluate.py --features proposed.csv --performance perf.csv \
        [--existing current.csv] [--objective min] [--good rel:0.05] \
        [--rho 0.1] [--pval 0.05]

CSV shape: first column is the instance id; any further non-numeric column is
ignored; every numeric column is a feature (or, in the performance file, an
algorithm). Ids must match across files.

The correlation check reproduces SIFTED's own filter, which is what decides
whether a candidate survives selection, so it runs Pearson on PRELIM-processed
values and takes its thresholds from the toolkit's own defaults. Against a built space,
read that space's model.opts.sifted and pass its rho and pval here instead of
trusting the defaults -- the run may have set them differently.

To audit a feature set rather than screen candidates against one, pass the set
itself as --features and leave --existing off. The degeneracy, correlation and
clustering reports then describe that set.

-------------------------------------------------------------------------
Instance Space Analysis (ISA) Toolkit
Copyright (c) 2026 Mario Andres Munoz Acosta and contributors
School of Computing and Information Systems
The University of Melbourne, Australia

SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
License: https://polyformproject.org/licenses/noncommercial/1.0.0/

You may use, modify, and redistribute this software for non-commercial
research and educational purposes only. Commercial use requires prior
written permission. See the LICENSE file for full terms.

Reference:
    Smith-Miles, K. & Munoz, M.A. (2023). Instance Space Analysis for
    Algorithm Testing. ACM Computing Surveys, 55(12), Article 255.
    https://doi.org/10.1145/3572895
"""

import argparse
import csv
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy import stats

NAN_THRESHOLD = 0.20  # opts.prelim.nanThreshold; INIT drops a column at or above it
IQR_MULTIPLIER = 5.0  # opts.prelim.iqrMultiplier

_LAMBDA_ZERO = 1e-10  # a Box-Cox lambda this small is the log transform
_MIN_ROWS = 3  # fewest observed rows a column needs to be transformed or correlated
_MIN_DISTINCT = 2  # fewest distinct values a column or label needs to correlate
_MIN_QUAD_ROWS = 10  # fewest rows, and distinct x values, for the quadratic fit
_MIN_QUAD_LEVELS = 3
_NEAR_CONSTANT = 0.90  # share of rows at the modal value that earns a WARN


def lab(name: str, w: int = 28) -> str:
    """Feature names run long; keep the tail, which is what distinguishes them."""
    return name if len(name) <= w else ".." + name[-(w - 2) :]


def read_csv(path: str) -> tuple[list[str], list[str], NDArray[np.float64]]:
    """Read an id column plus every numeric column, as an id list, names, matrix."""
    with Path(path).open() as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        sys.exit(f"{path}: empty")
    cols = list(rows[0])
    idcol = cols[0]
    num: list[str] = []
    for c in cols[1:]:
        try:
            [float(r[c] or "nan") for r in rows]
            num.append(c)
        except ValueError:
            pass
    ids = [r[idcol] for r in rows]
    x = np.array([[float(r[c] or "nan") for c in num] for r in rows], float)
    return ids, num, x


def align(
    ids_a: Sequence[str],
    ids_b: Sequence[str],
    x_b: NDArray[np.float64],
    what: str,
) -> NDArray[np.float64]:
    """Reorder the rows of x_b, keyed by ids_b, into the order of ids_a."""
    pos = {i: k for k, i in enumerate(ids_b)}
    missing = [i for i in ids_a if i not in pos]
    if missing:
        sys.exit(f"{what}: {len(missing)} ids absent, e.g. {missing[:3]}")
    return np.asarray(x_b[[pos[i] for i in ids_a]], dtype=float)


def boxcox(x: NDArray[np.float64], lam: float) -> NDArray[np.float64]:
    """PRELIM's Box-Cox application (isacompat.boxcoxApply)."""
    out = np.log(x) if abs(lam) < _LAMBDA_ZERO else np.expm1(lam * np.log(x)) / lam
    return np.asarray(out, dtype=float)


def boxcox_fit(x: NDArray[np.float64]) -> float:
    """Maximum-likelihood lambda on positive values; 1 for a constant sample."""
    if np.allclose(x, x[0]):
        return 1.0
    return float(stats.boxcox_normmax(x, method="mle"))


def zscore(x: NDArray[np.float64]) -> NDArray[np.float64]:
    """Standardise with the sample standard deviation, as PRELIM does."""
    mu, sigma = np.mean(x), np.std(x, ddof=1)
    return np.asarray((x - mu) / (sigma if sigma else 1.0), dtype=float)


def prelim_features(x: NDArray[np.float64]) -> NDArray[np.float64]:
    """Process feature columns the way PRELIM does.

    PRELIM's feature path: bound at median +/- 5*IQR, shift each column to a
    minimum of 1, Box-Cox, z-score. NaNs are left in place throughout.
    """
    med = np.nanmedian(x, 0)
    q75, q25 = np.nanpercentile(x, 75, axis=0), np.nanpercentile(x, 25, axis=0)
    iqr = q75 - q25
    x = np.clip(x, med - IQR_MULTIPLIER * iqr, med + IQR_MULTIPLIER * iqr)
    x = x - np.nanmin(x, 0) + 1.0
    out = np.full_like(x, np.nan)
    for j in range(x.shape[1]):
        ok = ~np.isnan(x[:, j])
        if ok.sum() < _MIN_ROWS:
            continue
        out[ok, j] = zscore(boxcox(x[ok, j], boxcox_fit(x[ok, j])))
    return out


def prelim_performance(y_perf: NDArray[np.float64]) -> NDArray[np.float64]:
    """Process performance columns the way PRELIM does.

    PRELIM's performance path. The shift is one global minimum over the whole
    matrix, not one per algorithm, and it adds eps rather than 1.
    """
    y_perf = (y_perf - np.nanmin(y_perf)) + np.finfo(float).eps
    out = np.full_like(y_perf, np.nan)
    for k in range(y_perf.shape[1]):
        ok = ~np.isnan(y_perf[:, k])
        if ok.sum() < _MIN_ROWS:
            continue
        out[ok, k] = zscore(boxcox(y_perf[ok, k], boxcox_fit(y_perf[ok, k])))
    return out


def quad_fit(x: NDArray[np.float64], y: NDArray[np.float64]) -> tuple[float, float]:
    """Multiple R for y ~ 1 + x + x^2 and its F-test p, pairwise-complete.

    Pearson scores a U-shape near zero, so a feature that predicts performance
    on both sides of a turning point fails the correlation check at nothing.
    This is what finds those before they are discarded in silence.
    """
    ok = ~np.isnan(x) & ~np.isnan(y)
    if ok.sum() < _MIN_QUAD_ROWS or len(np.unique(x[ok])) < _MIN_QUAD_LEVELS:
        return 0.0, 1.0
    xs, ys = x[ok], y[ok]
    if np.allclose(ys, ys[0]):
        return 0.0, 1.0
    design = np.column_stack([np.ones(xs.size), xs, xs**2])
    beta = np.linalg.lstsq(design, ys, rcond=None)[0]
    sstot = float(((ys - ys.mean()) ** 2).sum())
    if sstot <= 0:
        return 0.0, 1.0
    r2 = max(0.0, 1.0 - float(((ys - design @ beta) ** 2).sum()) / sstot)
    dfres = xs.size - 3
    if dfres <= 0:
        return 0.0, 1.0
    if r2 >= 1.0:
        return 1.0, 0.0
    f = (r2 / 2.0) / ((1.0 - r2) / dfres)
    return float(np.sqrt(r2)), float(stats.f.sf(f, 2, dfres))


def pearson_pairwise(
    x: NDArray[np.float64],
    y: NDArray[np.float64],
) -> tuple[float, float]:
    """corr(x, y, 'rows', 'pairwise'): drop a row if either value is missing."""
    ok = ~np.isnan(x) & ~np.isnan(y)
    if (
        ok.sum() < _MIN_ROWS
        or np.allclose(x[ok], x[ok][0])
        or np.allclose(y[ok], y[ok][0])
    ):
        return 0.0, 1.0
    r, p = stats.pearsonr(x[ok], y[ok])
    return (0.0, 1.0) if np.isnan(r) else (float(r), float(p))


def main() -> None:  # noqa: PLR0912
    """Run every check and print the report."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--features", required=True, help="candidate features CSV")
    ap.add_argument("--performance", required=True)
    ap.add_argument(
        "--existing",
        help="feature set already in use (enables the overlap columns)",
    )
    ap.add_argument(
        "--objective",
        choices=["min", "max"],
        default="min",
        help="min = opts.perf.MaxPerf false, the toolkit default",
    )
    ap.add_argument(
        "--good",
        default="rel:0.05",
        help="abs:<threshold> or rel:<epsilon>; mirrors AbsPerf + epsilon",
    )
    ap.add_argument(
        "--rho",
        type=float,
        default=0.10,
        help="opts.sifted.rho: minimum |correlation| to keep",
    )
    ap.add_argument(
        "--pval",
        type=float,
        default=0.05,
        help="opts.sifted.pval: significance level for that correlation",
    )
    ap.add_argument(
        "--K",
        type=int,
        default=10,
        help="opts.sifted.K: how many clusters SIFTED forms "
        "(default 10, from utils/ISAdefaults.m)",
    )
    a = ap.parse_args()

    ids, names, feat_raw = read_csv(a.features)
    pids, algos, perf_raw = read_csv(a.performance)
    perf_raw = align(ids, pids, perf_raw, a.performance)
    n, m = feat_raw.shape
    n_algos = len(algos)

    # ---- the good/bad rule, and the target SIFTED would correlate against ----
    try:
        kind, val_text = a.good.split(":")
        val = float(val_text)
        assert kind in ("abs", "rel")
    except (ValueError, AssertionError):
        sys.exit("--good must be abs:<threshold> or rel:<epsilon>")

    filler = np.inf if a.objective == "min" else -np.inf
    y_aux = np.where(np.isnan(perf_raw), filler, perf_raw)
    y_best = (y_aux.min(1) if a.objective == "min" else y_aux.max(1))[:, None]
    denom = np.where(y_best == 0, np.finfo(float).eps, y_best)
    y_nz = np.where(perf_raw == 0, np.finfo(float).eps, perf_raw)
    deficit = (y_nz / denom - 1.0) if a.objective == "min" else (1.0 - y_nz / denom)

    if kind == "abs":
        y_perf = perf_raw.copy()
        y_bin = (y_aux <= val) if a.objective == "min" else (y_aux >= val)
    else:
        y_perf = deficit
        y_bin = deficit <= val
    y_bin = y_bin.astype(float)

    print(f"n = {n} instances, {m} candidate features, {n_algos} algorithms")
    print(
        f"good rule: {kind} {val:g}, objective {a.objective}; "
        f"the correlation target is the "
        f"{'raw performance' if kind == 'abs' else 'deficit from the best'}",
    )
    print(
        "  good rate: "
        + ", ".join(f"{c} {y_bin[:, k].mean():.0%}" for k, c in enumerate(algos)),
    )
    dead = (y_bin.sum(1) == 0).mean()
    print(
        f"  instances no algorithm is good on: {dead:.0%}"
        + ("  <- these are ties in the deficit target" if kind == "rel" else "")
        + "\n",
    )

    # The build removes an algorithm that is good on no instance before PRELIM
    # normalises anything (InstanceSpace.m:636-645), so it is not in the matrix SIFTED
    # correlates against. Leave it out here too.
    never = y_bin.sum(0) == 0
    if never.all():
        sys.exit(
            "no algorithm is good on any instance under this rule, so the build "
            "would have no portfolio left",
        )
    if never.any():
        gone = [c for c, k in zip(algos, never) if k]
        print(
            "  never good, so the build removes them before anything else runs: "
            + ", ".join(gone)
            + ". They are left out below.\n",
        )
        y_perf, y_bin = y_perf[:, ~never], y_bin[:, ~never]
        algos = [c for c, k in zip(algos, never) if not k]
        n_algos = len(algos)

    # ---- reject 1: degeneracy -----------------------------------------------
    print("== reject 1: degeneracy ==")
    print(
        f"{'feature':28s} {'distinct':>8s} {'mode':>6s} {'missing':>8s} "
        f"{'minority':>9s}",
    )
    alive = np.ones(m, bool)
    for j, f in enumerate(names):
        v = feat_raw[:, j]
        obs = v[~np.isnan(v)]
        k = len(np.unique(obs))
        _, counts = np.unique(obs, return_counts=True)
        mode = counts.max() / len(obs) if len(obs) else 1.0
        nan = np.isnan(v).mean()
        iqr0 = len(obs) and np.percentile(obs, 75) == np.percentile(obs, 25)
        note, minority = "", "        -"
        if k <= 1:
            alive[j] = False
            note = "  FAIL constant"
        elif nan >= NAN_THRESHOLD:
            alive[j] = False
            note = f"  FAIL missing >= {NAN_THRESHOLD:.0%}, INIT would drop it"
        elif iqr0:
            alive[j] = False
            minority = f"{len(obs) - counts.max():9d}"
            note = "  FAIL zero IQR, PRELIM's bounding collapses it to a constant"
        elif mode >= _NEAR_CONSTANT:
            minority = f"{len(obs) - counts.max():9d}"
            note = (
                f"  WARN near-constant, {len(obs) - counts.max()} rows carry the signal"
            )
        print(f"{lab(f):28s} {k:8d} {mode:6.2f} {nan:8.2f} {minority}{note}")
    print(
        f"  {alive.sum()}/{m} pass. A WARN passes through to the next check; all it "
        "costs is\n  sample size, since the correlation rests on the minority rows "
        "alone.",
    )
    print(
        "  A zero IQR is different and is fatal as the pipeline is configured: the "
        "middle 50%\n  of the column is a single value, so median +/- 5*IQR has zero "
        "width and PRELIM's\n  bounding step (core/PRELIM.m:225) clips every "
        "row to "
        "that value. SIFTED then sees a\n  constant column, scores it 0 "
        "(core/SIFTED.m:117) and drops it. Any binary column\n  whose minority is "
        "under a quarter of the rows lands here. To keep such a feature,\n  set "
        "opts.bound.flag false, or redefine it so the middle 50% is not one value.\n",
    )

    # ---- reject 2: SIFTED's correlation filter ------------------------------
    feat = prelim_features(feat_raw)
    y_proc = prelim_performance(y_perf)
    r_mat = np.zeros((m, n_algos))
    p_mat = np.ones((m, n_algos))
    for j in range(m):
        for k in range(n_algos):
            r_mat[j, k], p_mat[j, k] = pearson_pairwise(feat[:, j], y_proc[:, k])

    print(
        f"== reject 2: Pearson on processed values, |rho| >= {a.rho} at p <= {a.pval} "
        "==",
    )
    print(
        "this is core/SIFTED.m:115 -- same target, same transform, same thresholds, "
        "no correction",
    )
    print(
        f"{'feature':28s} "
        + " ".join(f"{c[:9]:>10s}" for c in algos)
        + f" {'spread':>7s}",
    )
    passes = np.zeros(m, bool)
    for j, f in enumerate(names):
        hits = [
            abs(r_mat[j, k]) >= a.rho and p_mat[j, k] <= a.pval for k in range(n_algos)
        ]
        passes[j] = alive[j] and any(hits)
        cells = [f"{r_mat[j,k]:+9.3f}{'*' if hits[k] else ' '}" for k in range(n_algos)]
        print(
            f"{lab(f):28s} "
            + " ".join(cells)
            + f" {r_mat[j].max() - r_mat[j].min():7.3f}"
            + ("" if passes[j] else "  drop"),
        )
    print(
        f"  {passes.sum()}/{m} pass. '*' marks an algorithm this feature clears the "
        "floor on;\n  spread is the range across the portfolio, high meaning the "
        "feature separates\n  algorithms rather than only tracking difficulty.",
    )
    print(
        "  The floor is sufficient, not necessary: SIFTED also keeps whichever feature "
        "is\n  most correlated with each algorithm whatever its value "
        "(core/SIFTED.m:121), so in a\n  real run a candidate below the floor can "
        "still "
        "survive by being the best there is.",
    )

    curved = []
    for j in range(m):
        if passes[j] or not alive[j]:
            continue
        for k in range(n_algos):
            if abs(r_mat[j, k]) >= a.rho and p_mat[j, k] <= a.pval:
                continue
            rq, pq = quad_fit(feat[:, j], y_proc[:, k])
            if rq >= a.rho and pq <= a.pval and rq > abs(r_mat[j, k]):
                curved.append((j, k, rq))
    if curved:
        print(
            "\n  NON-MONOTONE SIGNAL in candidates the check above drops. A quadratic "
            "reaches\n  the same floor where Pearson does not, so the relationship is "
            "real and runs\n  both ways around a turning point:",
        )
        for j, k, rq in sorted(curved, key=lambda t: -t[2]):
            print(
                f"    {lab(names[j]):28s} {lab(algos[k], 20):20s} "
                f"Pearson {r_mat[j, k]:+.3f}, quadratic R {rq:.3f}",
            )
        print(
            "  The pipeline is Pearson throughout, so it will not find these as they\n"
            "  stand. Redefine the feature so the relationship runs one way -- the "
            "distance\n  from the turning point rather than the raw value -- and run "
            "it "
            "again.",
        )
    print()
    if passes.sum() == 0:
        sys.exit("nothing survives the correlation check")

    survivors = [names[j] for j in np.where(passes)[0]]

    # ---- rank ---------------------------------------------------------------
    rb_mat = np.zeros((m, n_algos))
    for j in range(m):
        for k in range(n_algos):
            ok = ~np.isnan(feat_raw[:, j])
            if (
                ok.sum() < _MIN_ROWS
                or len(np.unique(y_bin[ok, k])) < _MIN_DISTINCT
                or len(np.unique(feat_raw[ok, j])) < _MIN_DISTINCT
            ):
                continue
            r = stats.pointbiserialr(y_bin[ok, k], feat_raw[ok, j])[0]
            rb_mat[j, k] = 0.0 if np.isnan(r) else r

    exist: NDArray[np.float64] | None = None
    overlap: dict[str, tuple[float, str]] | None = None
    enames: list[str] = []
    if a.existing:
        eids, enames, exist_raw = read_csv(a.existing)
        exist = prelim_features(align(ids, eids, exist_raw, a.existing))
        overlap = {
            f: max(
                (
                    abs(pearson_pairwise(feat[:, names.index(f)], exist[:, k])[0]),
                    enames[k],
                )
                for k in range(exist.shape[1])
            )
            for f in survivors
        }

    print("== rank: point-biserial against the binary label, per algorithm ==")
    print("what PYTHIA predicts and TRACE draws footprints from, so it can order the")
    print("survivors differently from the check above")
    print(
        f"{'feature':28s} "
        + " ".join(f"{c[:9]:>10s}" for c in algos)
        + f" {'spread':>7s}",
    )
    for j in sorted(
        (names.index(f) for f in survivors),
        key=lambda j: -np.abs(rb_mat[j]).max(),
    ):
        print(
            f"{lab(names[j]):28s} "
            + " ".join(f"{rb_mat[j, k]:+10.3f}" for k in range(n_algos))
            + f" {rb_mat[j].max() - rb_mat[j].min():7.3f}",
        )
    print()

    if overlap is not None:
        print(f"== rank: overlap with {len(enames)} existing features ==")
        print(
            "whether the set already measures this. A duplicate does no harm, since "
            "SIFTED\nclusters the two and keeps one, but it buys nothing either",
        )
        print(f"{'feature':28s} {'max|r|':>7s}  {'nearest existing':28s}")
        for f in sorted(survivors, key=lambda f: -overlap[f][0]):
            r_near, c = overlap[f]
            print(f"{lab(f):28s} {r_near:7.3f}  {lab(c):28s}")
        print()

    pool = [(names[j], feat[:, j], "candidate") for j in np.where(passes)[0]]
    if exist is not None:
        pool += [(enames[k], exist[:, k], "existing") for k in range(exist.shape[1])]
    print("== rank: the clusters SIFTED would form ==")
    if len(pool) <= max(a.K, 3):
        print(
            f"  {len(pool)} column{'s' if len(pool) != 1 else ''} is not more than "
            f"K = {a.K}, so SIFTED skips the clustering\n  and keeps all of them "
            "(core/SIFTED.m:140).\n",
        )
    else:

        def tag(t: str) -> str:
            return "" if t == "candidate" or not a.existing else " (existing)"

        try:
            from sklearn.cluster import KMeans
        except ImportError:
            print("  (scikit-learn absent; skipping the clustering)\n")
        else:
            mat = np.vstack([c[1] for c in pool])
            mat = np.where(np.isnan(mat), np.nanmean(mat, 1, keepdims=True), mat)
            mat = mat - mat.mean(1, keepdims=True)
            nrm = np.linalg.norm(mat, axis=1, keepdims=True)
            mat = mat / np.where(nrm == 0, 1.0, nrm)
            grp = KMeans(a.K, n_init=10, random_state=0).fit_predict(mat)
            print(
                "  k-means on correlation distance into opts.sifted.K clusters, the way "
                "SIFTED\n  groups features (core/SIFTED.m:176). It then keeps "
                "exactly "
                "one column from\n  each (core/SIFTED.m:200-204), so columns "
                "sharing a "
                "cluster are competing and\n  all but one of them will be dropped "
                "whatever you add to the set.",
            )
            for g in range(a.K):
                members = [pool[i] for i in np.where(grp == g)[0]]
                if not members:
                    continue
                print(
                    f"    cluster {g + 1}: "
                    + ", ".join(lab(nm, 24) + tag(kind) for nm, _, kind in members),
                )
            print(
                f"  {len(pool)} columns in {len(set(grp))} clusters, so the set comes "
                "out "
                f"of selection\n  at {len(set(grp))} columns however many go in.\n",
            )

    if exist is not None:
        try:
            from sklearn.impute import SimpleImputer  # type: ignore[import-untyped]
            from sklearn.linear_model import LogisticRegression
            from sklearn.model_selection import StratifiedKFold, cross_val_score
            from sklearn.pipeline import make_pipeline  # type: ignore[import-untyped]
        except ImportError:
            print("  (scikit-learn absent; skipping incremental value)\n")
            return
        x_existing = exist
        x_combined = np.hstack([exist, feat[:, [names.index(f) for f in survivors]]])
        cv = StratifiedKFold(5, shuffle=True, random_state=0)

        def mk() -> Any:  # noqa: ANN401
            return make_pipeline(SimpleImputer(), LogisticRegression(max_iter=5000))

        print(
            "== rank: incremental value, 5-fold CV AUC, existing vs existing + "
            "survivors ==",
        )
        print(
            f"{'algorithm':16s} {'existing':>9s} {'+survivors':>11s} {'gain':>8s} "
            f"{'headroom':>9s}",
        )
        for k, c in enumerate(algos):
            y = y_bin[:, k]
            n_minority = int(min(y.sum(), len(y) - y.sum()))
            if len(np.unique(y)) < _MIN_DISTINCT:
                print(f"{lab(c, 16):16s}  single class, skipped")
                continue
            if n_minority < cv.get_n_splits():
                print(
                    f"{lab(c, 16):16s}  {n_minority} in the smaller class, too few for "
                    f"{cv.get_n_splits()}-fold, skipped",
                )
                continue
            s1 = cross_val_score(mk(), x_existing, y, cv=cv, scoring="roc_auc").mean()
            s2 = cross_val_score(mk(), x_combined, y, cv=cv, scoring="roc_auc").mean()
            print(
                f"{lab(c, 16):16s} {s1:9.3f} {s2:11.3f} {s2 - s1:+8.3f} "
                f"{1.0 - s1:9.3f}",
            )
        print(
            "  read the gain against the headroom beside it: a small gain on an "
            "already-high\n  baseline is a ceiling effect, not evidence the candidates "
            "add nothing\n",
        )

    print("== rank: ordered ==")
    print(
        "on |rho| against the performance values, because that is what decides whether",
    )
    print(
        "SIFTED keeps the feature. 'on' names the algorithm each figure belongs to; "
        "the",
    )
    print(
        "two can be different algorithms, and either column can order the set its own",
    )
    print(
        "way. Discount a candidate whose overlap above is near 1 however high its rho "
        "--",
    )
    print("the order below does not apply that discount for you.")
    print(
        f"{'feature':28s} {'|rho| Y':>8s} {'on':11s} {'|r| Ybin':>9s} {'on':11s} "
        f"{'spread':>7s}",
    )
    for j in sorted(
        (names.index(f) for f in survivors),
        key=lambda j: -np.abs(r_mat[j]).max(),
    ):
        ky, kb = int(np.argmax(np.abs(r_mat[j]))), int(np.argmax(np.abs(rb_mat[j])))
        print(
            f"{lab(names[j]):28s} {np.abs(r_mat[j]).max():8.3f} "
            f"{lab(algos[ky], 11):11s} "
            f"{np.abs(rb_mat[j]).max():9.3f} {lab(algos[kb], 11):11s} "
            f"{r_mat[j].max() - r_mat[j].min():7.3f}",
        )

    # ---- yield --------------------------------------------------------------
    print("\n== yield ==")
    print(
        f"  {len(survivors)} of {m} candidates evaluated in this run survive both "
        "checks.",
    )
    print("  Report that against every candidate ever evaluated, not just this run's.")


if __name__ == "__main__":
    main()
