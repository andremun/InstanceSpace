r"""Performance-only diagnostics for ISA good-performance thresholds.

Answers one question: if you set ``opts.perf`` to these values, what labels will
PRELIM actually produce, and are they any good?

Reads a standard ISA metadata CSV and takes its ``algo_*`` columns, plus a
``Source`` column when one is present.

Usage as a library::

    from threshold_sweep import threshold_sweep
    r = threshold_sweep("metadata.csv", max_perf=False, abs_perf=True,
                        epsilon=[1, 10, 60, 600, 3599])
    print(r.data)      # threshold-independent degeneracy report -- read first
    print(r.sweep)     # one row per candidate rung

Usage from the shell::

    python threshold_sweep.py metadata.csv --minimise --absolute \
        --epsilon 1 10 60 600 3599

The binarisation replicates PRELIM's semantics exactly, including its NaN
substitution (-inf maximising, +inf minimising) and its zero-best-score eps
substitution. This is a faithful copy rather than an improvement: the job is to
*predict* what PRELIM will produce, quirks included, so a cleaner
reimplementation would be a worse tool. Verify with::

    [~, ~, out] = PRELIM(X, Y, opts);      % in MATLAB
    assert (binarise(Y, ...) == out.Ybin).all()

Companion script for the performance-thresholds skill. Replicates the good/bad
labelling of the Instance Space Analysis (ISA) Toolkit's core/PRELIM.m.

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

from __future__ import annotations

import argparse
import re
import warnings
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike, NDArray

__all__ = ["threshold_sweep", "binarise", "best_per_instance", "SweepResult"]

EPS = np.finfo(float).eps  # MATLAB's `eps`, for the zero-best substitution

# Good-rates outside [_RATE_LOW, _RATE_HIGH] count as saturated or near-degenerate.
_RATE_LOW = 0.05
_RATE_HIGH = 0.95
# PRELIM's cut for too many zero best-scores, and the share of values a point
# mass needs before it is treated as a sentinel.
_ZERO_BEST_FRAC = 0.05
_POINT_MASS_FRAC = 0.05
# An algorithm uniquely good on less than this share of instances is a passenger.
_PASSENGER_FRAC = 0.01
# Fewest points for a rank correlation, and fewest algorithms for a screen score.
_MIN_RANK_POINTS = 3
_MIN_SCREEN_ALGOS = 3
_MIN_PAIR_ALGOS = 2


# ---------------------------------------------------------------------------
# Core labelling -- a deliberate line-for-line port of core/PRELIM.m lines 85-127
# ---------------------------------------------------------------------------
def best_per_instance(y: NDArray[np.float64], max_perf: bool) -> NDArray[np.float64]:
    """Per-instance best score, matching PRELIM's NaN handling."""
    aux = np.where(np.isnan(y), -np.inf if max_perf else np.inf, y)
    return np.asarray(aux.max(axis=1) if max_perf else aux.min(axis=1), dtype=float)


def binarise(
    y: NDArray[np.float64],
    max_perf: bool,
    abs_perf: bool,
    epsilon: float,
) -> NDArray[np.bool_]:
    """Good/bad labels, replicating core/PRELIM.m exactly.

    Mirrors PRELIM rather than improving on it. Note the relative branch divides
    the NaN-substituted matrix by the *eps-substituted* best, and that PRELIM's
    separate ``Y(Y==0) = eps`` applies only to its returned Y, not to Ybin.
    """
    if max_perf:
        aux = np.where(np.isnan(y), -np.inf, y)
        ybest = aux.max(axis=1)
        if abs_perf:
            return np.asarray(aux >= epsilon, dtype=bool)
        ybest = np.where(ybest == 0, EPS, ybest)
        with np.errstate(invalid="ignore"):
            return np.asarray((1.0 - aux / ybest[:, None]) <= epsilon, dtype=bool)
    else:
        aux = np.where(np.isnan(y), np.inf, y)
        ybest = aux.min(axis=1)
        if abs_perf:
            return np.asarray(aux <= epsilon, dtype=bool)
        ybest = np.where(ybest == 0, EPS, ybest)
        with np.errstate(invalid="ignore"):
            return np.asarray((aux / ybest[:, None] - 1.0) <= epsilon, dtype=bool)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _read_algo_columns(
    csv_path: str,
) -> tuple[NDArray[np.float64], list[str], NDArray[np.str_] | None]:
    """Pull the algo_* columns (and any Source column) out of an ISA CSV."""
    df = pd.read_csv(csv_path)
    algo_cols = [c for c in df.columns if str(c).lower().startswith("algo_")]
    if not algo_cols:
        raise ValueError(f'No columns starting with "algo_" found in {csv_path}.')
    block = df[algo_cols]
    if not all(pd.api.types.is_numeric_dtype(block[c]) for c in algo_cols):
        raise ValueError(
            "At least one algo_ column is non-numeric. Check the CSV for "
            '"NA", blank rows, or Excel error codes such as #REF! or #DIV/0!.',
        )
    y = block.to_numpy(dtype=float)
    labels = [re.sub(r"^algo_", "", str(c), flags=re.IGNORECASE) for c in algo_cols]
    # The Source column, when present, names each instance's benchmark class or
    # generator. It is the population the degeneracy guard checks against -- a
    # threshold that merely re-derives it has sorted instances, not algorithms.
    src = next((c for c in df.columns if str(c).lower() == "source"), None)
    groups = df[src].astype(str).to_numpy() if src is not None else None
    return y, labels, groups


def _group_index(
    groups: NDArray[np.str_] | None,
) -> tuple[NDArray[np.intp] | None, int, bool]:
    """Map group labels to an ordered index.

    Benchmark classes usually carry a trailing number whose numeric order is the
    meaningful one -- alphabetical would put ``set10`` before ``set2``.
    ``ordered`` is False when the labels are not numbered, or when two families
    share a numbering (``setA02`` and ``setB02`` both parse to 2); groupRho is
    then meaningless and returned as NaN, leaving groupSpread as the guard.
    """
    if groups is None:
        return None, 0, False
    uniq, idx = np.unique(groups, return_inverse=True)
    nums: list[int | None] = []
    for u in uniq:
        m = re.search(r"(\d+)\s*$", str(u))
        nums.append(int(m.group(1)) if m else None)
    if all(n is not None for n in nums) and len(set(nums)) == len(nums):
        order = np.argsort(np.asarray(nums))
        remap = np.empty(len(order), dtype=int)
        remap[order] = np.arange(len(order))
        return remap[idx], len(uniq), True
    return idx, len(uniq), False


def _spearman(x: ArrayLike, y: ArrayLike) -> float:
    if len(np.atleast_1d(x)) < _MIN_RANK_POINTS:
        return np.nan
    xs, ys = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    r = pd.Series(xs).corr(pd.Series(ys), method="spearman")
    return float(r) if r is not None and np.isfinite(r) else np.nan


def _mean_pairwise_jaccard(lab: NDArray[np.bool_]) -> float:
    """Mean overlap between algorithms' good-sets.

    Lower means the portfolio is more complementary at this threshold. Pairs
    whose union is empty carry no information and are skipped.
    """
    n = lab.shape[1]
    if n < _MIN_PAIR_ALGOS:
        return np.nan
    vals: list[float] = []
    for a in range(n - 1):
        for b in range(a + 1, n):
            union = np.count_nonzero(lab[:, a] | lab[:, b])
            if union:
                vals.append(np.count_nonzero(lab[:, a] & lab[:, b]) / union)
    return float(np.mean(vals)) if vals else np.nan


def _plateau_stats(
    grid: NDArray[np.float64],
    g_rate: NDArray[np.float64],
    e: float,
    g_at_e: float,
    tol: float,
) -> tuple[float, float]:
    """Local slope and plateau width around a candidate rung.

    Width is measured on the dense curve: the widest contiguous interval
    containing ``e`` over which the good-rate stays within ``tol`` of its value
    at ``e``. A wide plateau is the evidence that the exact value is not doing
    the work.
    """
    slope: float = np.nan
    width: float = np.nan
    if e < grid[0] or e > grid[-1]:
        return slope, width
    i = int(np.argmin(np.abs(grid - e)))
    if 0 < i < len(grid) - 1:
        slope = float((g_rate[i + 1] - g_rate[i - 1]) / (grid[i + 1] - grid[i - 1]))
    # A flat curve at a saturated good-rate is degeneracy, not robustness: the
    # labels have stopped changing, so of course epsilon can move freely. Report
    # no width rather than one that would be misread as evidence of stability.
    if g_at_e < _RATE_LOW or g_at_e > _RATE_HIGH:
        return slope, np.nan
    ok = np.abs(g_rate - g_at_e) <= tol
    if not ok[i]:
        return slope, 0.0
    lo = hi = i
    while lo > 0 and ok[lo - 1]:
        lo -= 1
    while hi < len(grid) - 1 and ok[hi + 1]:
        hi += 1
    return slope, float(grid[hi] - grid[lo])


def _ladder_125(lo: float, hi: float) -> NDArray[np.float64]:
    """Fallback 1-2-5 sweep grid spanning the data range.

    Exists so the script runs without a caller-supplied ladder. It is a sweep
    convenience, not a proposal -- a real ladder comes from the metric's units
    and the domain's conventions (see the skill's "Generate the candidate
    ladder", §3).
    """
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= 0:
        return np.array([0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0], dtype=float)
    lo = max(lo, hi / 1e4)
    decades = range(int(np.floor(np.log10(lo))), int(np.ceil(np.log10(hi))) + 1)
    e = np.unique(np.array([m * 10.0**d for d in decades for m in (1, 2, 5)]))
    e = e[(e >= lo / 2) & (e <= hi * 1.5)]
    return e if e.size else np.linspace(lo, hi, 7)


# ---------------------------------------------------------------------------
@dataclass
class SweepResult:
    """Container for a sweep. See the skill for how to read each column."""

    data: dict[str, Any]  # threshold-independent degeneracy report
    per_algo: pd.DataFrame  # per-algorithm min/max/median/range
    sweep: pd.DataFrame  # one row per candidate rung
    curve: pd.DataFrame  # dense good-rate curve
    opts: dict[str, Any]
    labels: list[NDArray[np.bool_]] = field(default_factory=list)  # Ybin per rung
    unique_by_algo: list[NDArray[np.float64]] = field(
        default_factory=list,
    )  # per-algo unique share


def threshold_sweep(  # noqa: PLR0912
    csv_path: str,
    *,
    max_perf: bool = False,
    abs_perf: bool = False,
    epsilon: Sequence[float] | None = None,
    plateau_tol: float = 0.05,
    grid_points: int = 200,
) -> SweepResult:
    """Sweep candidate good-performance thresholds over an ISA metadata CSV."""
    y, algo_labels, groups = _read_algo_columns(csv_path)
    ninst, nalgos = y.shape
    finite = y[np.isfinite(y)]
    ybest = best_per_instance(y, max_perf)

    uniq_vals = np.unique(finite)
    resolution = float(np.diff(uniq_vals).min()) if uniq_vals.size > 1 else np.nan
    if finite.size:
        vals, counts = np.unique(finite, return_counts=True)
        pm_val, pm_frac = float(vals[counts.argmax()]), counts.max() / finite.size
    else:
        pm_val, pm_frac = np.nan, 0.0

    # A sentinel (timeout, ceiling) needs a rung strictly inside it, and how far
    # inside is not free: any genuine value between the rung and the sentinel is
    # relabelled as a censored run. These are the nearest genuine values on each
    # side of the point mass, so the caller can see the interval a rung may
    # safely occupy instead of guessing at cutoff - 1.
    below = uniq_vals[uniq_vals < pm_val] if uniq_vals.size else np.array([])
    above = uniq_vals[uniq_vals > pm_val] if uniq_vals.size else np.array([])
    pm_below = float(below.max()) if below.size else np.nan
    pm_above = float(above.min()) if above.size else np.nan

    # Ties matter whenever epsilon = 0 is a live candidate: at that rung every
    # tied-best algorithm is labelled good, so the tie rate is the entire
    # deviation from a clean winner-take-all split. Mirrors the multiple-best
    # statistic PRELIM reports.
    tie_frac = float(np.mean((ybest[:, None] == y).sum(axis=1) > 1))

    data: dict[str, Any] = {
        "nInstances": ninst,
        "nAlgos": nalgos,
        "algoLabels": algo_labels,
        "globalMin": float(finite.min()) if finite.size else np.nan,
        "globalMax": float(finite.max()) if finite.size else np.nan,
        "nanFrac": float(np.mean(np.isnan(y))),
        "ybestZeroFrac": float(np.mean(ybest == 0)),
        "tieFrac": tie_frac,
        "hasNegative": bool((finite < 0).any()),
        "crossesZero": bool((finite < 0).any() and (finite > 0).any()),
        "pointMassValue": pm_val,
        "pointMassFrac": float(pm_frac),
        "pointMassNeighbourBelow": pm_below,
        "pointMassNeighbourAbove": pm_above,
        "resolution": resolution,
        "nDistinct": int(uniq_vals.size),
    }

    # When the point mass sits at an extreme it is almost always a sentinel, and
    # the rung has to clear it without swallowing genuine runs. Name the safe
    # interval explicitly -- the failure it prevents is silent.
    # A sentinel that must be labelled bad has to sit on the bad side of the
    # rung, which pins it to the max under minimising and to the min under
    # maximising; those are the only two arrangements a rung can resolve. A mass
    # on the *good* side cannot be cleared by any epsilon, so say that instead of
    # returning nothing.
    if pm_frac >= _POINT_MASS_FRAC and np.isfinite(pm_val):
        at_max, at_min = pm_val == data["globalMax"], pm_val == data["globalMin"]
        if at_max and not max_perf:
            data["safeRungInterval"] = (pm_below, pm_val)  # minimising: (gap, sentinel)
        elif at_min and max_perf:
            data["safeRungInterval"] = (pm_val, pm_above)  # maximising: (sentinel, gap)
        elif at_max or at_min:
            data["safeRungInterval"] = "goodSide"
        else:
            data["safeRungInterval"] = None

    # Relative mode is meaningless on a scale that hits or crosses zero. PRELIM
    # substitutes eps and warns at the same 5% cut; raise it here too so the
    # caller sees it before proposing a mode.
    if not abs_perf and data["ybestZeroFrac"] > _ZERO_BEST_FRAC:
        warnings.warn(
            f'{100 * data["ybestZeroFrac"]:.1f}% of instances have a '
            "best-algorithm score of exactly zero. Relative mode is degenerate "
            "on this data -- use abs_perf=True.",
            RuntimeWarning,
            stacklevel=2,
        )
    if not abs_perf and data["crossesZero"]:
        warnings.warn(
            "Performance values cross zero, so the relative ratio has no "
            "consistent meaning. Use abs_perf=True.",
            RuntimeWarning,
            stacklevel=2,
        )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        per_algo = pd.DataFrame(
            {
                "algorithm": algo_labels,
                "min": np.nanmin(y, axis=0),
                "max": np.nanmax(y, axis=0),
                "median": np.nanmedian(y, axis=0),
                "nanFrac": np.isnan(y).mean(axis=0),
            },
        )
    per_algo["range"] = per_algo["max"] - per_algo["min"]

    # Dense good-rate curve, computed independently of the reported rungs so
    # plateau width and local slope are not limited by a coarse ladder.
    if abs_perf:
        lo, hi = data["globalMin"], data["globalMax"]
        if not np.isfinite(lo) or not np.isfinite(hi) or lo == hi:
            lo, hi = 0.0, 1.0
        pad = 0.02 * (hi - lo)
        grid = np.linspace(lo - pad, hi + pad, grid_points)
    else:
        grid = np.linspace(0.0, 1.0, grid_points)
    g_rate = np.array([binarise(y, max_perf, abs_perf, e).mean() for e in grid])
    curve = pd.DataFrame({"epsilon": grid, "goodRate": g_rate})

    # Candidate rungs. The default is an explicit fallback for sweeping, NOT a
    # proposal -- the ladder is supposed to come from the domain (skill §3).
    ladder: ArrayLike = (
        epsilon
        if epsilon is not None
        else (
            _ladder_125(data["globalMin"], data["globalMax"])
            if abs_perf
            else np.array([0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5])
        )
    )
    eps = np.sort(np.atleast_1d(np.asarray(ladder, dtype=float)))

    # Absolute rungs are in the metric's own units; relative rungs are
    # fractions of the per-instance best. Passing a relative ladder
    # (0.01, 0.05, ...) while asking for absolute mode is the common slip, and
    # it fails silently -- every rung below the data's minimum labels
    # everything bad, every rung above its maximum labels everything good.
    if abs_perf and np.isfinite(data["globalMin"]):
        stray = eps[(eps < data["globalMin"]) | (eps > data["globalMax"])]
        if stray.size:
            warnings.warn(
                f"{stray.size} of {eps.size} absolute rungs fall outside the "
                f'data range [{data["globalMin"]:.6g}, {data["globalMax"]:.6g}]: '
                f"{np.array2string(stray, precision=6)}. Absolute rungs are in "
                "the metric's own units. If you meant fractions of the "
                "per-instance best, use abs_perf=False.",
                RuntimeWarning,
                stacklevel=2,
            )

    gidx, ngroups, ordered = _group_index(groups)
    rows: list[dict[str, Any]] = []
    all_labels: list[NDArray[np.bool_]] = []
    uniq_by_algo: list[NDArray[np.float64]] = []

    for e in eps:
        lab = binarise(y, max_perf, abs_perf, e)
        all_labels.append(lab)
        good_rate = float(lab.mean())
        algo_rate = lab.mean(axis=0)
        ngood = lab.sum(axis=1)
        frac_all_good = float(np.mean(ngood == nalgos))
        frac_all_bad = float(np.mean(ngood == 0))
        jac = _mean_pairwise_jaccard(lab)

        # Screening components. Balance peaks at a 50% good-rate;
        # informativeness is the share of instances that discriminate between at
        # least two algorithms; complementarity rewards non-overlapping sets.
        balance = 1 - abs(2 * good_rate - 1)
        inform = 1 - frac_all_good - frac_all_bad
        compl = 0.0 if np.isnan(jac) else 1 - jac
        screen = float(np.cbrt(max(balance, 0) * max(inform, 0) * max(compl, 0)))
        # With two algorithms fracAllBad is structurally zero and meanJaccard
        # equals fracAllGood, so the score collapses to a monotone function of
        # epsilon peaking at 0 -- it measures disagreement, not discrimination.
        if nalgos < _MIN_SCREEN_ALGOS:
            screen = np.nan

        # Per-algorithm unique contribution: instances this algorithm is the
        # only one to get right. An algorithm with a decent good-rate but
        # near-zero unique contribution is a passenger.
        uniq = np.array(
            [
                float(np.mean(lab[:, a] & ~np.delete(lab, a, axis=1).any(axis=1)))
                for a in range(nalgos)
            ],
        )
        uniq_by_algo.append(uniq)

        # Degeneracy guard: does this rung sort instance classes rather than
        # algorithms? Rejection criterion only -- never an objective.
        g_spread = g_rho = np.nan
        if ngroups > 1:
            gr = np.array([lab[gidx == g].mean() for g in range(ngroups)])
            g_spread = float(gr.max() - gr.min())
            if ordered:
                g_rho = _spearman(np.arange(ngroups), gr)

        slope, width = _plateau_stats(grid, g_rate, float(e), good_rate, plateau_tol)
        rows.append(
            {
                "epsilon": e,
                "goodRate": good_rate,
                "algoRateMin": float(algo_rate.min()),
                "algoRateMax": float(algo_rate.max()),
                "nDegenerateAlgos": int(np.sum((algo_rate == 0) | (algo_rate == 1))),
                "fracAllGood": frac_all_good,
                "fracAllBad": frac_all_bad,
                "meanJaccard": jac,
                "uniqueBestFrac": float(np.mean(ngood == 1)),
                "groupSpread": g_spread,
                "groupRho": g_rho,
                # Exact 0/1 misses an algorithm good on 91% or 5% of instances,
                # which is still a badly unbalanced target for PYTHIA.
                # Informational, not disqualifying.
                "nNearDegenerate": int(
                    np.sum((algo_rate <= _RATE_LOW) | (algo_rate >= _RATE_HIGH)),
                ),
                "nPassengers": int(np.sum(uniq < _PASSENGER_FRAC)),
                "dGoodRate": slope,
                "plateauWidth": width,
                "screenScore": screen,
            },
        )

    return SweepResult(
        data=data,
        per_algo=per_algo,
        sweep=pd.DataFrame(rows),
        curve=curve,
        opts={
            "max_perf": max_perf,
            "abs_perf": abs_perf,
            "plateau_tol": plateau_tol,
            "grid_points": grid_points,
        },
        labels=all_labels,
        unique_by_algo=uniq_by_algo,
    )


# ---------------------------------------------------------------------------
def _fmt(x: float) -> str:
    return f"{x:.3f}"


def _main() -> None:
    ap = argparse.ArgumentParser(description=(__doc__ or "").split("\n")[0])
    ap.add_argument("csv")
    d = ap.add_mutually_exclusive_group()
    d.add_argument("--minimise", action="store_true", help="lower is better (default)")
    d.add_argument("--maximise", action="store_true", help="higher is better")
    m = ap.add_mutually_exclusive_group()
    m.add_argument(
        "--absolute",
        action="store_true",
        help="fixed threshold in metric units",
    )
    m.add_argument(
        "--relative",
        action="store_true",
        help="fraction of the per-instance best (default)",
    )
    ap.add_argument(
        "--epsilon",
        type=float,
        nargs="+",
        default=None,
        help="candidate thresholds -- in the metric's own units when "
        "--absolute, as fractions of the per-instance best "
        "when --relative (default: a fallback sweep grid)",
    )
    a = ap.parse_args()

    r = threshold_sweep(
        a.csv,
        max_perf=a.maximise,
        abs_perf=a.absolute,
        epsilon=a.epsilon,
    )

    print("== data (threshold-independent; read this first) ==")
    for k, v in r.data.items():
        shown = v
        if k == "algoLabels":
            shown = ", ".join(v)
        if k == "safeRungInterval":
            if v is None:
                continue
            if v == "goodSide":
                print(
                    f"  {k:<24}point mass sits on the GOOD side of any value: "
                    "if it is a sentinel, no epsilon labels it bad",
                )
                continue
            lo, hi = v
            if not np.isfinite(lo) or not np.isfinite(hi):
                shown = "no genuine values on the safe side of the sentinel"
            else:
                shown = (
                    f"a threshold strictly between {lo:g} and {hi:g} labels the "
                    f"sentinel bad and every genuine run good"
                )
        print(f"  {k:<24}{shown}")
    print("\n== per algorithm ==")
    print(r.per_algo.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    # Two blocks rather than one, both keyed by epsilon so rows still line up
    # across them. Every r.sweep column the skill's references name as a
    # diagnostic has to be visible here, and all of them on one line is 200
    # characters.
    # epsilon gets %g, not the fixed 3 decimals the diagnostics get: candidate
    # values span 0.0001 to 3600 depending on the metric and the mode, and
    # rounding would print 0 and 0.0001 identically -- two candidates that mean
    # quite different things.
    eps_fmt: Mapping[str | int, Callable[..., Any]] = {"epsilon": lambda x: f"{x:g}"}
    print("\n== sweep: labels ==")
    print(
        r.sweep[
            [
                "epsilon",
                "goodRate",
                "algoRateMin",
                "algoRateMax",
                "fracAllGood",
                "fracAllBad",
                "nDegenerateAlgos",
                "nNearDegenerate",
            ]
        ].to_string(index=False, float_format=_fmt, formatters=eps_fmt),
    )
    print("\n== sweep: selection and shape ==")
    print(
        r.sweep[
            [
                "epsilon",
                "meanJaccard",
                "uniqueBestFrac",
                "nPassengers",
                "dGoodRate",
                "plateauWidth",
                "groupSpread",
                "groupRho",
                "screenScore",
            ]
        ].to_string(index=False, float_format=_fmt, formatters=eps_fmt),
    )
    print("\n== unique contribution per algorithm ==")
    for e, u in zip(r.sweep["epsilon"], r.unique_by_algo):
        joined = "  ".join(f"{n}={x:.3f}" for n, x in zip(r.data["algoLabels"], u))
        print(f"  eps={e:<12g}{joined}")


if __name__ == "__main__":
    _main()
