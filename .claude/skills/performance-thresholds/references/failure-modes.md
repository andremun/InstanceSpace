# Failure modes

Five ways a threshold decision goes wrong, each with the diagnostic that
catches it and what to do about it. `SKILL.md` lists the symptoms. Read
the entry that matches what you are seeing; if nothing has gone wrong,
you need none of this.

Field names here come from `scripts/threshold_sweep.m`;
`sweep-columns.md` says what each one is.

## No candidate value is acceptable

The good-rate jumps between two neighbouring candidates — 5% labels
almost nothing good, 10% almost everything. Before reaching for
something in between, check whether a round value sits in the gap, 7.5%
here. If one does, it belonged on the shortlist: add it and test it.

Otherwise the cliff has a cause, and it will be one of these:

- a point mass (check `pointMassFrac`) — usually a timeout, so read
  *Censored metrics* below;
- bimodal instance difficulty — genuinely two populations, and worth
  saying so;
- a degenerate portfolio, one algorithm dominating (check
  `algoRateMax`);
- a metric that needs a transform, such as log-runtime, before any
  threshold means much.

Once you know which, an intermediate value may still be the best you can
do. If so, recommend it: say it came from the sweep and not from the
domain, and give the range over which the labels stay the same, so
nobody takes it for a domain choice.

## Relative mode on a scale that hits or crosses zero

PRELIM substitutes `eps` for zero best-scores and warns above 5%
(`ISA:PRELIM:manyZeroBest`), but the resulting ratios are meaningless.
Check `ybestZeroFrac`, `hasNegative` and `crossesZero` in `r.data` before
proposing relative mode. Classification error where the best classifier
reaches exactly 0.0 is the common case.

## Censored metrics

Two traps, both fatal and neither raising an error. The cutoff itself is
usually the right absolute threshold and needs no justifying — "solved
within the budget" means something on its own. But:

**1. It cannot be written as `epsilon = cutoff`.** PRELIM's absolute test
is `Y <= epsilon` when minimising, and censored runs are recorded at
*exactly* the cutoff — so `epsilon = B` on a budget of `B` marks every
censored run good, every algorithm good on every instance, and the label
carries no information. Encode it strictly inside the sentinel. The same
holds when maximising against a ceiling, where the test is
`Y >= epsilon` and the value must sit strictly above the capped one. This
produces no error, just a meaningless space.

**Backing off from the sentinel is not free.** `cutoff - 1` is the
reflex, and it is only correct when nothing genuine sits in the interval
it skips: any real run recorded there is silently relabelled as censored.
On a 3600 s budget, `epsilon = 3599` misfiles a solver that returned at
3599.4 s as a timeout — and soft time limits, coarse clock polls, and
solvers that check their budget between nodes all produce values in that
last second. So **never propose a backed-off value without checking the
gap**: `safeRungInterval` in `r.data` states the range a value may occupy
and `pointMassNeighbourBelow` the value it must clear, both printed by
the sweep. Put it anywhere inside that interval and report the interval,
since the labels are identical across all of it. If the interval is
empty, or narrower than the metric's recorded precision, no absolute
threshold separates the sentinel from real data; say so and ask for a
run-status flag.

**2. Relative mode inverts on censored data.** Where every algorithm hit
the limit, the per-instance best *is* the sentinel, so every algorithm
sits at deficit 0 and is labelled good. Relative mode then asserts
universal success precisely on the instances where the portfolio
universally failed, and no `epsilon` fixes it. Penalising (PAR10) does
not rescue it either — scaling the sentinel preserves the tie among
censored runs, so the inversion survives. A censored metric with any
meaningful all-fail mass goes to absolute mode.

The sentinel also **misorders, not just censors**: a run that crashed at
40 s is recorded identically to one still progressing at the limit. If
the sentinel absorbs failures as well as timeouts, the bad label mixes
two mechanisms — ask for a run-status flag, and report the split.

## A relative threshold quantises on a coarse integer metric

PRELIM's relative test admits about `epsilon * best` units of slack,
truncated to whole units, so on a count-valued metric — features
selected, errors made, messages sent, resources allocated — the slack a
given `epsilon` actually grants varies with the instance. That much is
the point of relative mode, and it is not the defect: a tolerance that
scales with the per-instance best is usually what you meant by a
percentage, and on a metric whose values run into the thousands the
truncation is invisible.

The defect is the collapse at the bottom of the range. Once
`floor(epsilon * best)` reaches zero, those instances get a tie test —
only an algorithm matching the best exactly is good — while the largest
instances are still being given tens of units. Nobody chose the tie test,
and the reported percentage describes neither group. Worse, instance
scale is usually family-determined, so the tolerance ends up assigned
**by benchmark family**: the label partitions instance classes instead of
algorithms, which is what `groupSpread` catches, and no `epsilon` fixes
it.

So the question is not whether the metric is an integer, it is whether
the smallest instances still get slack. Check `resolution` against the
*smallest* per-instance best, not the mean: if `epsilon * best` stays
comfortably above `resolution` there, relative mode is sound and nothing
below applies. If it does not, you want *"within k units of the best"*
instead. That is **not expressible with a scalar `epsilon`** — pass the
excess matrix and use absolute mode:

```matlab
D = Y - min(Y, [], 2);        % or max(...) when maximising
opts.perf.AbsPerf = true;
opts.perf.epsilon = k;
```

Write `D` into the `algo_` columns of `metadata.csv` before the build.

`k = 0` is the one case where the two agree: absolute-0 on `D` and
relative-0 on `Y` produce identical labels, so "matches the best found"
needs no transform.

Two things to know about `D`. It fixes the quantisation, but it is
portfolio-derived, so a space built on it cannot carry a hardness claim
however absolute the mode looks. And it leaves a large point mass at 0
with a long right tail, which PRELIM's Box-Cox step handles poorly —
check the normalised result before trusting anything built on it.

## Hardness mapping needs an external bound

No `epsilon` substitutes for one. Any threshold measured against the
*portfolio's own* best — relative mode, or absolute mode on
`Y - min(Y)` — makes the best algorithm good by construction, so all-bad
instances cannot exist. If the question is "where does the whole
portfolio fail", it needs a bound from outside the portfolio: a proven
optimum, a relaxation or combinatorial bound, a published best-known
value, an irreducible-error estimate, a human or reference baseline. Say
the question is unanswerable on the data rather than shipping a
portfolio-relative space described as a hardness map. Partial coverage is
often enough — a bound column plus a coverage flag lets you build the
space on the covered subset.
