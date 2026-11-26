# `threshold_sweep.m` output, field by field

Load this when you have sweep output in front of you. `SKILL.md` carries
the rules — which fields rule out a mode, what never to tune against —
and this is the glossary behind them. Nothing here decides anything; if a
field looks like it is telling you what to choose, the skill is where
that gets settled.

## `r.data` — the threshold-independent report

Computed before any threshold, and identical whether you ask for absolute
or relative — which is what lets it rule out a mode before you commit to
one. It does depend on `MaxPerf`, since `tieFrac` and `ybestZeroFrac` are
computed from the per-instance best. Read it first:

- per-algorithm min / max / median / NaN fraction, and the global range
  (feeds the candidates in §3);
- `tieFrac` — fraction of instances where more than one algorithm attains
  the best score. Only matters when `epsilon = 0` is in play, where it is
  the entire deviation from a clean winner-take-all split;
- `ybestZeroFrac` — fraction of instances whose best score is exactly
  zero. **Above 0.05 this kills relative mode**; PRELIM raises
  `ISA:PRELIM:manyZeroBest` at the same cut;
- `hasNegative`, `crossesZero` — either makes the relative ratio
  meaningless;
- `pointMassValue` / `pointMassFrac` — censoring and timeout detection;
- `pointMassNeighbourBelow` / `pointMassNeighbourAbove` and
  `safeRungInterval` — the nearest genuine values either side of that mass,
  and the range an absolute candidate may occupy without swallowing them.
  **Read this before proposing any value near a sentinel**; `cutoff - 1` is
  only
  safe when this interval says it is;
- `resolution`, `nDistinct` — the smallest gap between distinct values, a
  floor on how finely any threshold can meaningfully cut.
Three of those settle things on their own — `ybestZeroFrac`,
`hasNegative`/`crossesZero`, and `safeRungInterval` — and are repeated as
rules in `SKILL.md`, because they decide things and the skill has to
carry them whether or not this file is loaded.

## `r.sweep` — one row per candidate value

| Column | Read it for |
|---|---|
| `goodRate` | Overall balance. Near 0 or 1 is degenerate. |
| `algoRateMin`, `algoRateMax`, `nDegenerateAlgos` | An algorithm at 100% good contributes nothing. One at 0% is worse: the build removes any algorithm with no good instances before anything else runs (`InstanceSpace.m:636-645`), so it disappears from every later table and figure. |
| `nNearDegenerate` | Algorithms at or beyond 5%/95%. Informational, not disqualifying — a dominant algorithm is a finding — but it warns that PYTHIA's target is badly unbalanced for those algorithms. |
| `fracAllGood`, `fracAllBad` | Instances that discriminate between *no* algorithms. Note `fracAllBad` is structurally 0 in relative mode. |
| `meanJaccard` | Overlap between algorithms' good-sets. **Lower is better** — it means the portfolio is complementary at this threshold. |
| `uniqueBestFrac` | Instances where exactly one algorithm is good. Direct measure of selection value. |
| `groupSpread`, `groupRho` | **The degeneracy guard.** Spread of the good-rate across instance groups (a `Source` column, if present), and its rank correlation with group index — `NaN` when the labels carry no unique ordering. Read it against algorithm discrimination, per the rule below. Reject on this; never tune on it. |
| `nPassengers` | How many algorithms are the *only* good one on almost no instances. A passenger inflates `fracAllBad` while adding nothing to selection. The per-algorithm figures behind it are in `r.uniqueByAlgo`, not in this table. |
| `dGoodRate` | Local slope of the good-rate. Large magnitude = cliff. A central difference, so `NaN` at either end of the swept range. |
| `plateauWidth` | How far `epsilon` can move before `goodRate` shifts by more than 0.05. Wide = the exact value is not doing the work. `NaN` when the good-rate is already saturated — see below. |
| `screenScore` | See the warning below. `NaN` with fewer than three algorithms. |
## Why `screenScore` is `NaN` for a two-algorithm portfolio

With two algorithms `fracAllBad` is structurally zero and `meanJaccard`
equals `fracAllGood`, so the score collapses to a monotone function of
`epsilon` peaking at 0 — it measures disagreement rather than
discrimination. Read `uniqueBestFrac` and the per-algorithm rates
instead.

## Why `plateauWidth` is `NaN` at a saturated good-rate

It is returned as `NaN` whenever `goodRate` is below 0.05 or above 0.95.
A flat curve there is not robustness: if 97% of pairs are already
labelled good the labels have stopped changing, so of course the curve is
flat. Returning a width would invite reading degeneracy as evidence.

## The other outputs

- `r.perAlgo` — per-algorithm ranges, for generating candidates.
- `r.curve` — the dense good-rate curve behind `plateauWidth` and
  `dGoodRate`.
- `r.labels{k}` — the actual `Ybin` for candidate `k`, for checking
  anything this script does not compute.
- `r.uniqueByAlgo{k}` — for candidate `k`, each algorithm's unique
  contribution: the share of instances where it is the only good one.
  `r.sweep`'s `nPassengers` counts how many of these fall below 0.01.

## A note on `safeRungInterval`

The name is the script's and predates the skill dropping the word "rung"
from its prose. It is the interval a candidate threshold may occupy near
a sentinel without swallowing real data. Renaming it would break callers,
so it stays.
