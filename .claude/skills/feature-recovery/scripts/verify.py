"""Verify candidate feature definitions against the column values.

You write definitions.py; this checks whether each one actually reproduces its
column, and when it does not, says what kind of wrong it is.

    definitions.py must provide:
        load(path) -> {instance id: parsed instance}
        DEFINITIONS = {"feature_name": callable(instance) -> float, ...}

    load may optionally take an `ids` keyword; when --sample is in use it is
    given the ids actually being verified, so only those instances need parsing.

    python verify.py --features features.csv --defs definitions.py \
        [--instances PATH] [--id-col NAME] [--quiet] \
        [--columns A,B] [--sample N] [--sample-seed K]

While searching for a reading, narrow the run rather than widening it:
--columns to the one column under study, --sample to 20 or 30 rows. A
subsample refutes a candidate; only the full run confirms one.

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
"""

import argparse
import csv
import importlib.util
import inspect
import random
import sys
import time
from collections import deque
from collections.abc import Iterable, Sequence
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
from numpy.typing import NDArray

# ------------------------------------------------------------------ tolerance --


def decimals(s: str) -> int | None:
    """Decimal places a CSV cell was printed with, or None for scientific notation."""
    s = s.strip()
    return None if "e" in s.lower() else (len(s.split(".")[1]) if "." in s else 0)


def numeric(s: str) -> bool:
    """Whether a CSV cell parses as a number."""
    try:
        float(s)
        return True
    except ValueError:
        return False


# A column stored as float32 before it was printed carries about seven
# significant digits however many decimals the CSV shows, so a definition that
# is right still lands this far out.  Storage width cannot be read back
# reliably once a CSV has rounded to a fixed number of decimal *places*, so
# rather than guess it, every column is scored at both bounds and told apart by
# which one it meets.
SINGLE = 1e-6

# Thresholds the diagnosis uses to name a near miss.
_NOISE_REL = 1e-10  # last-bits disagreement between two exact computations
_EXACT_FIT = 1e-12  # correlation, slope and offset indistinguishable from exact
_UNIT_SLOPE = 1e-9  # a fitted slope this close to 1 is no rescaling at all
_NEAR_PROPORTIONAL = 0.999999
_SAME_ORDER = 0.9999
# Columns listed by name before the rest are summarised as "...".
_SHOW_UNDEFINED = 6
_MIN_SPLIT = 2  # an interval shorter than this has no interior point to bisect


def tolerances(strings: Sequence[str], is_int: bool) -> NDArray[np.float64]:
    """Return the rounding tolerance for each cell of a column.

    Per value, not per column: a CSV at fixed significant figures carries
    less absolute precision at larger magnitudes, so a column-wide bound
    wrongly rejects the big rows. Slack keeps an exact half-ulp tie passing.
    """
    if is_int:
        return np.full(len(strings), 1e-9)
    eps = float(np.finfo(float).eps)
    out: list[float] = []
    for s in strings:
        x = abs(float(s))
        d = decimals(s)
        rounded = 1e-9 * max(1.0, x) if d is None else 0.5 * 10.0**-d * (1 + 1e-6)
        # A CSV written at full double precision prints the exact double, so the
        # rounding bound collapses below one ulp.  Never demand better than float
        # reproduction: recomputing reorders operations and moves the last bits.
        out.append(max(rounded, 8 * eps * max(1.0, x)))
    return np.array(out)


def loose(y: NDArray[np.float64]) -> NDArray[np.float64]:
    """Tolerance at single precision, for telling rounding from a wrong reading."""
    return np.asarray(SINGLE * np.maximum(1.0, np.abs(y)), dtype=float)


# ---------------------------------------------------------------- diagnostics --


def partial_note(hit: int, okf: NDArray[np.bool_], n: int) -> str:
    """Describe a partial match at both tolerances.

    Report both bounds; the failures are diagnosed at the looser one, since
    rows that differ only by rounding are not the ones to go and characterise.
    """
    note = f"matches {hit}/{n}"
    if int(okf.sum()) > hit:
        note += (
            f"; {int(okf.sum())}/{n} to ~7 significant digits, so "
            f"{n - int(okf.sum())} are real disagreements and the rest "
            f"is rounding"
        )
    return note


def diagnose(  # noqa: PLR0911, PLR0912
    y: NDArray[np.float64],
    v: NDArray[np.float64],
    tol: NDArray[np.float64],
) -> tuple[str, str, NDArray[np.bool_] | None]:
    """Classify how a definition's values compare with the column.

    Name what kind of wrong a near miss is. Each verdict points somewhere
    specific, so none of them are 'try again'.
    """
    ok = np.abs(v - y) <= tol
    n, hit = len(y), int(ok.sum())
    if hit == n:
        return "EXACT", "", None
    if not np.all(np.isfinite(v)):
        return "ERROR", "definition returned non-finite values", None
    if len(np.unique(np.round(v, 12))) == 1 and len(np.unique(np.round(y, 12))) > 1:
        return "NO MATCH", "definition is constant; the column is not", ~ok

    # Same quantity, last bits apart: iterative solvers and long accumulations
    # cannot reproduce bit-for-bit, so this is a pass with a note rather than a
    # miss -- kept distinct from a genuine scale error.
    scale = np.maximum(1.0, np.abs(y))
    if np.max(np.abs(v - y) / scale) < _NOISE_REL:
        return (
            "EXACT~",
            f"within numerical noise (max |Δ| = {np.max(np.abs(v - y)):.2g})",
            None,
        )

    # Every row agrees to single precision.  A column printed with more decimals
    # than it was stored with cannot do better, and neither can a long
    # accumulation, so the reading is right and the residue is rounding.
    okf = np.abs(v - y) <= loose(y)
    if okf.all():
        return (
            "EXACT~f32",
            f"agrees to ~7 significant digits (max rel = "
            f"{np.max(np.abs(v - y) / scale):.1e}); the column is not stored "
            f"more precisely than this",
            None,
        )

    # A reading exact on most rows is missing a branch, not a transform, however
    # well it correlates -- a few far-out rows drag the fit below.  Count only
    # the non-zero rows: a wrong scale matches every zero for free.
    nz = np.abs(y) > tol
    if nz.any() and int((okf & nz).sum()) * 2 > int(nz.sum()):
        return "PARTIAL", partial_note(hit, okf, n), ~okf

    # affine: right quantity, wrong scale -- the fitted slope names the gap
    if np.std(v) > 0:
        a, b = np.polyfit(v, y, 1)
        r = np.corrcoef(v, y)[0, 1]
        if r > 1 - _EXACT_FIT:
            hint = f"scale by {a:.10g}" + (
                f", offset {b:+.6g}" if abs(b) > _EXACT_FIT else ""
            )
            if abs(a - 1) > _UNIT_SLOPE or abs(b) > _EXACT_FIT:
                return "AFFINE", f"exactly proportional -- {hint}", ~ok
        if r > _NEAR_PROPORTIONAL:
            return "AFFINE", f"near-proportional (r={r:.9f}), scale {a:.8g}", ~ok
        from scipy.stats import spearmanr

        rho = spearmanr(v, y).statistic
        if rho > _SAME_ORDER:
            return "MONOTONE", f"same ordering (rho={rho:.6f}); unknown transform", ~ok
    if hit:
        return "PARTIAL", partial_note(hit, okf, n), ~okf
    if okf.any():
        return "PARTIAL", f"matches {int(okf.sum())}/{n} to ~7 significant digits", ~okf
    return "NO MATCH", "no relationship", ~ok


def structured(
    bad: NDArray[np.bool_],
    ids: Sequence[str],
    y: NDArray[np.float64],
    v: NDArray[np.float64],
) -> str:
    """Summarise the failing rows of a column.

    Scattered failures are edge cases; a predictable set is a missing
    branch. Reported with the failing ids so the pattern can be characterised
    -- and weigh any claimed pattern against how common it is overall.
    """
    k = int(bad.sum())
    if k == 0:
        return ""
    ex = [ids[i] for i in np.where(bad)[0][:4]]
    d = np.abs(v[bad] - y[bad])
    return (
        f"      {k} failing, e.g. {ex}  |Δ| median {np.median(d):.3g} max {d.max():.3g}"
    )


# -------------------------------------------------------------------- sample --


def bisect_order(n: int) -> list[int]:
    """Positions 0..n-1 ordered ends-first, then repeatedly bisecting.

    Yields a nested sequence -- the first k of it is a sensible spread for any
    k -- so raising --sample extends the sample instead of replacing it, and a
    row that refuted a candidate at 20 is still there at 50.
    """
    if n <= 0:
        return []
    if n == 1:
        return [0]
    out: list[int] = []
    seen: set[int] = set()
    q = deque([(0, n - 1)])
    for p in (0, n - 1):
        out.append(p)
        seen.add(p)
    while q:
        lo, hi = q.popleft()
        if hi - lo < _MIN_SPLIT:
            continue
        mid = (lo + hi) // 2
        if mid not in seen:
            out.append(mid)
            seen.add(mid)
        q.append((lo, mid))
        q.append((mid, hi))
    return out


def value_lane(
    rows: Sequence[dict[str, str]],
    idx: Iterable[int],
    c: str,
) -> list[int] | None:
    """Rows of `idx` ordered ends-first by column c, or None if c is not numeric."""
    try:
        vals = {i: float(rows[i][c]) for i in idx}
    except (ValueError, KeyError):
        return None
    order = sorted(idx, key=lambda i: vals[i])
    return [order[p] for p in bisect_order(len(order))]


def choose_sample(
    rows: Sequence[dict[str, str]],
    cols: Sequence[str],
    n: int,
    seed: int,
    groups: Sequence[str] | None = None,
) -> tuple[list[int], int]:
    """Choose the rows a --sample run verifies.

    Pick n rows spread across the value range of each column under study,
    and across the sources if the CSV names one per row.

    Not a uniform draw.  What a subsample has to catch is a missing branch, and
    those sit at the ends of a column's range -- the degenerate instance, the
    one row exactly on a boundary, the largest instance -- which random
    sampling reaches with low probability precisely because they are rare.

    Two strata, in priority order.  The global extremes of each column go in
    first, since a boundary convention shows up there and nowhere else.  Then
    the picks are spread across the sources rather than left to cluster in one,
    because the other common missing branch is a whole generator handled
    differently.  Within a source, rows are again taken ends-first.

    Deterministic, so two candidate readings are compared on the same rows.
    """
    m = len(rows)
    if n >= m:
        return list(range(m)), 0

    picks: list[int] = []
    seen: set[int] = set()

    def take(i: int) -> None:
        if i not in seen and len(picks) < n:
            picks.append(i)
            seen.add(i)

    col_lanes = [lane for lane in (value_lane(rows, range(m), c) for c in cols) if lane]

    # The extremes first -- but capped, so that on a wide CSV they cannot eat
    # the whole budget and crowd out the sources.
    for lane in col_lanes:
        if len(picks) >= max(2, n // 2):
            break
        for i in lane[:2]:
            take(i)

    # One row per source before any source gets a second.
    src_lanes: list[list[int]] = []
    if groups:
        bysrc: dict[str, list[int]] = {}
        for i, g in enumerate(groups):
            bysrc.setdefault(g, []).append(i)
        for g in sorted(bysrc):
            src_lane = value_lane(rows, bysrc[g], cols[0]) if cols else None
            src_lanes.append(src_lane or bysrc[g])

    for lanes in (src_lanes, col_lanes):
        for tier in range(max((len(lane) for lane in lanes), default=0)):
            if len(picks) >= n:
                break
            for lane in lanes:
                if tier < len(lane):
                    take(lane[tier])
    spread = len(picks)

    rest = [i for i in range(m) if i not in seen]
    random.Random(seed).shuffle(rest)
    picks += rest[: max(0, n - len(picks))]
    return sorted(picks), spread


def call_load(
    mod: ModuleType,
    path: str,
    ids: set[str] | None,
) -> tuple[dict[str, Any], bool]:
    """Pass the sampled ids to load() when it accepts them.

    Parsing is often the expensive half, and with --sample most instances are
    not wanted.  Optional, so definitions.py written against the two-argument
    contract keeps working.
    """
    params: Any
    try:
        params = inspect.signature(mod.load).parameters
    except (TypeError, ValueError):
        params = {}
    takes = "ids" in params or any(
        p.kind is inspect.Parameter.VAR_KEYWORD for p in params.values()
    )
    if takes and ids is not None:
        loaded: dict[str, Any] = mod.load(path, ids=ids)
        return loaded, True
    loaded = mod.load(path)
    return loaded, False


# ---------------------------------------------------------------------- main --


def main() -> None:  # noqa: PLR0912
    """Verify every requested column and print one verdict per column."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--features", required=True)
    ap.add_argument("--defs", required=True)
    ap.add_argument("--instances", default=None)
    ap.add_argument(
        "--id-col",
        default=None,
        help="instance id column (default: the first column)",
    )
    ap.add_argument("--quiet", action="store_true", help="only show what is not EXACT")
    ap.add_argument(
        "--columns",
        default=None,
        help="comma-separated columns to verify (default: all of them); "
        "narrow to the column under study while searching",
    )
    ap.add_argument(
        "--sample",
        type=int,
        default=None,
        metavar="N",
        help="verify on N rows drawn from the column extremes and "
        "spread across the source column, for the search loop; "
        "20-30 is the working size. Refutes a reading, never "
        "confirms one",
    )
    ap.add_argument(
        "--sample-seed",
        type=int,
        default=0,
        help="seed for the non-extreme part of the sample (default: 0)",
    )
    a = ap.parse_args()

    spec = importlib.util.spec_from_file_location("defs", a.defs)
    if spec is None or spec.loader is None:
        sys.exit(f"cannot import {a.defs}")
    mod = importlib.util.module_from_spec(spec)
    # Its own directory goes on the path first, so definitions.py can split
    # helpers -- an expensive kernel, a wrapper round a compiled one -- into
    # sibling modules instead of being one file by necessity.
    sys.path.insert(0, str(Path(a.defs).resolve().parent))
    spec.loader.exec_module(mod)
    for need in ("load", "DEFINITIONS"):
        if not hasattr(mod, need):
            sys.exit(f"{a.defs} must define {need}")

    with Path(a.features).open() as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        sys.exit(f"{a.features} has no data rows")
    idcol = a.id_col or next(iter(rows[0]))
    if idcol not in rows[0]:
        sys.exit(f"no column {idcol!r}; columns are {list(rows[0])[:8]}")
    cols = [c for c in rows[0] if c != idcol and c.lower() not in ("source", "sources")]
    if a.columns:
        want = [c.strip() for c in a.columns.split(",") if c.strip()]
        missing = [c for c in want if c not in cols]
        if missing:
            sys.exit(f"no such column(s): {missing}; columns are {cols[:8]}")
        cols = want

    # Whether each column varies over the *whole* CSV, recorded before any
    # subsetting: a column that is flat only because the sample is small must
    # not be reported as CONSTANT, which in this skill means unidentifiable.
    varies = {c: len({r[c] for r in rows}) > 1 for c in cols}
    full_rows = rows

    # Sample against the columns actually being verified, so the spread is over
    # the ranges under study rather than over columns nobody asked about.
    # The source column is excluded from scoring but is the best stratum there
    # is: a missing branch is often one generator handled differently.
    gcol = next((c for c in rows[0] if c.lower() in ("source", "sources")), None)
    picked = None
    if a.sample is not None and a.sample < len(rows):
        groups = [r[gcol] for r in rows] if gcol else None
        idx, spread = choose_sample(rows, cols, a.sample, a.sample_seed, groups)
        rows = [rows[i] for i in idx]
        picked = (len(idx), spread)

    nsrc = len({r[gcol] for r in full_rows}) if gcol else 0
    ids_all = [r[idcol] for r in rows]
    objs, narrowed = call_load(
        mod,
        a.instances or str(Path(a.features).resolve().parent),
        set(ids_all) if picked else None,
    )

    ids = ids_all
    keep = [i for i, x in enumerate(ids) if x in objs]
    if picked:
        n, spread = picked
        srcs = len({r[gcol] for r in rows}) if gcol else 0
        covered = f", covering {srcs} of {nsrc} sources" if gcol else ""
        print(
            f"sample    : {n} rows, {spread} taken from the ends and quantiles "
            f"of {len(cols)} column(s){covered}"
            + (
                ""
                if narrowed
                else "; load() parsed every instance -- give it an "
                "`ids` keyword to parse only these"
            ),
        )
        print("            a sample refutes a reading; only the full run confirms one")
    else:
        print(
            f"alignment : {len(keep)}/{len(rows)} rows have instances; "
            f"{len(objs) - len(keep)} instances unused",
        )
    if not keep:
        sys.exit("nothing aligned -- check the id convention first")
    if len(keep) < len(rows):
        print(
            f"            verification denominator is {len(keep)}, not {len(rows)}; "
            "unmatched rows can refute a definition but never confirm one",
        )
    kids = [ids[i] for i in keep]

    undefined = [c for c in cols if c not in mod.DEFINITIONS]
    print(f"defined   : {len(cols) - len(undefined)}/{len(cols)} columns\n")

    tally: dict[str, int] = {}
    clock: dict[str, float] = {}

    def fail(col: str, msg: str) -> None:
        """One column's problem must not end the run: report it and carry on."""
        print(f"{col:<52} ERROR       {msg}")
        tally["ERROR"] = tally.get("ERROR", 0) + 1

    for c in cols:
        if c not in mod.DEFINITIONS:
            continue
        s = [rows[i][c].strip() for i in keep]
        try:
            y = np.array([float(x) for x in s])
        except ValueError:
            bad_cells = [x for x in s if not numeric(x)]
            fail(c, f"{len(bad_cells)} non-numeric cell(s), e.g. {bad_cells[:3]}")
            continue
        if not np.all(np.isfinite(y)):
            # NaN/inf parse cleanly and would then silently poison every
            # comparison, so stop here rather than report a meaningless verdict.
            fail(c, f"{int((~np.isfinite(y)).sum())} non-finite cell(s) in the column")
            continue
        is_int = bool(np.all(y == np.round(y))) and all(decimals(x) == 0 for x in s)
        tol = tolerances(s, is_int)
        t0 = time.perf_counter()
        try:
            v = np.array([float(mod.DEFINITIONS[c](objs[i])) for i in kids])
        except Exception as e:
            fail(c, f"{type(e).__name__}: {e}")
            continue
        clock[c] = time.perf_counter() - t0
        if len(np.unique(np.round(y, 12))) == 1:
            if picked and varies.get(c):
                print(
                    f"{c:<52} FLAT        constant in this sample but not in the "
                    f"column; raise --sample to score it",
                )
                tally["FLAT"] = tally.get("FLAT", 0) + 1
                continue
            print(f"{c:<52} CONSTANT    column never varies; any constant matches it")
            tally["CONSTANT"] = tally.get("CONSTANT", 0) + 1
            continue
        verdict, note, bad = diagnose(y, v, tol)
        tally[verdict] = tally.get(verdict, 0) + 1
        if verdict in ("EXACT", "EXACT~", "EXACT~f32") and a.quiet:
            continue
        secs = f"  [{clock[c]:.1f}s]" if clock[c] >= 1.0 else ""
        print(f"{c:<52} {verdict:<11} {note}{secs}")
        if bad is not None and verdict != "NO MATCH":
            print(structured(bad, kids, y, v))

    print("\n" + "  ".join(f"{k}: {v}" for k, v in sorted(tally.items())))
    if undefined:
        print(
            f"no definition supplied for {len(undefined)}: {undefined[:6]}"
            + (" ..." if len(undefined) > _SHOW_UNDEFINED else ""),
        )

    # Which column to go and make cheaper -- a guess about where the time went
    # is usually wrong, and the answer decides whether --sample is worth it.
    total = sum(clock.values())
    if total >= 1.0:
        # Only the columns worth acting on: a list ending in "0.0s" buries the
        # one that matters.
        slow = [
            kv
            for kv in sorted(clock.items(), key=lambda kv: -kv[1])
            if kv[1] >= 0.05 * total
        ][:3]
        per = f"{total / max(1, len(keep)) * 1e3:.1f} ms/instance"
        print(
            f"time      : {total:.1f}s over {len(keep)} instances ({per}); slowest "
            + ", ".join(f"{c} {t:.1f}s" for c, t in slow),
        )


if __name__ == "__main__":
    main()
