---
name: performance-thresholds
description: >
  Propose and justify the "good performance" threshold that ISA's PRELIM
  stage uses to binarise algorithm performance — the absolute-vs-relative
  choice and the threshold value. Use this skill whenever the user asks
  what threshold to use, what counts as good performance, how to set
  `opts.perf.epsilon` / `AbsPerf` / `MaxPerf`, whether a threshold is
  discriminating or degenerate, or asks to justify or sensitivity-test an
  existing one. Trigger on "performance threshold", "what counts as
  good", "epsilon", "AbsPerf", "absolute or relative", "binarise
  performance", "good/bad label", or a request to review the threshold
  behind an existing instance space. Covers the method, the diagnostics,
  and `scripts/threshold_sweep.m`. For the surrounding methodology see
  the `instance-space-analysis` skill.
---

# Performance thresholds for Instance Space Analysis

Choose the mode and the threshold value PRELIM uses to turn performance
numbers into good/bad labels, and write down the argument that defends
them.

If a space has already been built and the threshold is only one suspect
among several, `instance-space-improvement` works out whether it is the
live problem before you get here. For the surrounding methodology see
`instance-space-analysis`.

Two reference files, loaded when needed rather than on every run:
`references/sweep-columns.md` for what each sweep output means, and
`references/failure-modes.md` for the five ways this goes wrong.

---

## 1. The decision

Three knobs, in `opts.perf`, applied by `core/PRELIM.m`:

| Knob | Meaning |
|---|---|
| `MaxPerf` | Direction. `true` = maximise (accuracy, profit), `false` = minimise (error, runtime, cost). Almost always given by the metric, not chosen. |
| `AbsPerf` | **The structural choice.** `true` = fixed threshold, `false` = relative to the per-instance best. |
| `epsilon` | The threshold value. Metric units when absolute; a *fraction* of the best when relative (`0.05` = within 5%). |

Every recommendation sets all three knobs explicitly rather than
relying on the defaults (*Pipeline defaults* in
`instance-space-interpretation`). An existing space's threshold is
read from its own `model.opts.perf`, never assumed from a default.

### Absolute vs relative is not a tuning choice

The two modes answer different questions, and the difference is
structural, not a matter of degree:

- **Relative** answers *"who wins here?"* The per-instance best algorithm
  has deficit 0, and `0 <= epsilon` always holds, so **the best algorithm
  is always labelled good**. There are never all-bad instances. You
  therefore cannot identify regions that defeat the entire portfolio. It
  also needs a metric that stays strictly positive: a scale that reaches
  or crosses zero makes the ratio meaningless.
- **Absolute** answers *"is this instance solved?"* It permits both
  all-good (trivially easy) and all-bad (genuinely hard) instances, which
  is what lets you locate regions that defeat the whole portfolio.

Decide this from the **intended use of the space** before touching
`epsilon`:

| Intended use | Mode |
|---|---|
| Algorithm selection / portfolio construction | Relative |
| Identifying hard regions, benchmark adequacy, "where does everything fail" | Absolute, **against a bound from outside the portfolio** (see below) |
| Claiming an instance is solved, or met a stated standard | Absolute |
| Comparing algorithms whose absolute scale is not meaningful or varies per instance | Relative |
| **The achievable maximum varies per instance** and no per-instance bound is available | **Relative** (see below) |

Settle the mode first. An answer that goes straight to `epsilon` has
skipped the decision that matters more.

**Absolute mode is necessary for a hardness claim, but not sufficient.**
What makes all-bad instances possible is a threshold that does not depend
on the portfolio. Absolute mode on a portfolio-derived quantity does not
qualify: on `Y - min(Y)` the best algorithm sits at a deficit of zero, so
it is good at any `epsilon >= 0` and nothing is ever all-bad — the same
construction as relative mode, wearing absolute mode's clothes. Note that
`references/failure-modes.md` recommends that very transform as the fix
for a *different* problem, a relative threshold quantising on a coarse
integer metric, so the two uses have to be kept apart.

### "Normalised" does not mean "comparable across instances"

A metric can be bounded and unitless and still have an achievable maximum
that differs from one instance to the next: a utilisation or fill
fraction, where each instance's structure caps what *any* solution can
reach, or an accuracy limited by irreducible error.

A fixed cut on such a metric labels instances by their own structure
rather than by algorithm quality. And because benchmark suites are
usually designed to vary exactly that structure, the label ends up
tracking the suite's design axis — the space re-derives the benchmark's
class index and says nothing about algorithms.

**Rule: a normalised metric with an instance-dependent maximum and no
per-instance bound goes to relative mode.** Where per-instance bounds or
proven optima can be had, better still: convert the metric to a
gap-to-bound and use absolute mode on that. Check
`groupSpread`/`groupRho` (§4) before committing to absolute on a metric
like this.

---

## 2. What you need to know

Ask for what is missing. Do not infer these from the performance matrix.

**Domain**
- What the problem is; what the algorithms are.
- What the literature in this subfield reports as "good", "solved" or
  "acceptable" — whatever that field's own success criterion is called.
- **Intended use of the instance space** (per the table above).

**Metric semantics**
- Direction, and units.
- Is zero attainable? Does the scale cross zero or go negative? Either
  kills relative mode.
- Bounded (accuracy, AUC, F1) or unbounded (runtime, objective value)?
- Is it a gap/ratio (already normalised) or a raw value in problem units?
- **Censoring**: is there a timeout or cutoff, and what is it? Is the
  metric penalised (PAR10)?
- **Measurement noise**: repeated runs per (algorithm, instance) pair, or
  a stated precision. Without one you cannot tell whether a candidate
  separates signal from rounding.

**Portfolio**
- Which algorithms are deliberately weak baselines vs. current SOTA.
- Best-known or proven-optimal values per instance, if available — these
  enable gap-based absolute thresholds, which are usually far more
  defensible than raw-value ones.

---

## 3. Propose the candidates

**The diagnostics filter candidates; the domain justification selects
them.** That is the rule the rest of this skill rests on, and it is why
the candidates come first: read the sweep before you have named any, and
you will anchor on wherever `screenScore` happens to peak — which is
never a reason for a choice.

### Where the candidates come from

The metric's units and range, not a formula.

**The two modes take values on different scales, and mixing them up is
the easiest mistake to make here.** Relative candidates are always
fractions of the per-instance best — 0, 1%, 2%, 5%, 10%, whatever the
metric measures, because the units cancel. Absolute candidates are always
in the metric's own units: seconds for a runtime, an error rate for a
classifier, MB for memory, whole units for anything counted.
`epsilon = 0.05`
is a *5% error rate* on a metric bounded in [0, 1]; on a runtime in
seconds it is fifty milliseconds, and on a count it is below the smallest
possible value, so every algorithm is labelled bad everywhere. If a
proposed absolute value falls outside the data's own range, it is almost
always a relative candidate applied to the wrong mode. The sweep warns
about this.

**Span the data.** Hour-scale values are absurd when every run finishes
in under 2 s; 1% is absurd when every gap is 20-80%. The sweep reports
per-algorithm min/max/range for exactly this (`r.perAlgo`). Looking at the
*range* before you have named any candidates is safe: a range is not the
good-rate curve, so it cannot tell you where `screenScore` peaks.

**Use the metric's own notion of round.** On a runtime in seconds, 3600
is round because it is an hour; 1000 is not, though it is rounder in
base ten. Time goes in human units, memory in powers of two, MIP gaps in
0.01/0.1/1, frame budgets in 16.7/33 ms. No rounding function recovers
this — it is the part of the job that needs domain knowledge.

Illustrative only, not a lookup table:

| Metric | Plausible candidates |
|---|---|
| Runtime, seconds | 1, 5, 10, 30, 60, 300, 600, 1800, 3600 |
| Optimality gap | 0%, 0.01%, 0.1%, 1%, 2%, 5%, 10% |
| Misclassification error | 1%, 5%, 10%, 20%, 25% |
| Bounded score (AUC) | 0.5 (chance), 0.7, 0.8, 0.9, 0.95 |
| Relative (`AbsPerf=false`) | 0%, 1%, 2%, 5%, 10%, 20% |

`epsilon = 0` in relative mode is often the right answer rather than an
edge case: it means "strictly the best", ties aside. It is the roundest
value available and the one that lets least instance structure leak into
the label, and it belongs on every relative shortlist. Check `tieFrac`
before recommending it, and check that the metric's measurement noise is
genuinely smaller than the winning margins.

### Which candidate to prefer

**One you could have named before seeing the results.** That is what
makes a value a commitment rather than a tuned parameter, and it is the
whole of the criterion. The strongest are values that mean something
independent of your data — the timeout itself, 0% gap, chance level, the
majority-class rate, a real-time or SLA budget — followed by whatever
this subfield already reports. Roundness is only a proxy for the same
thing, so a hard domain constraint beats it: a 33 ms frame budget or a
250 ms SLA is ugly and wins anyway.

Write the mode down with its argument from intended use, a shortlist of
2-4 values, and why each one. Nobody but you reads this, so it can be
long and it can hedge. Write it before running the sweep, so that what
you compare against is what you actually thought beforehand.

---

## 4. Check them against the data

`scripts/threshold_sweep.m` is a MATLAB function that also runs in GNU
Octave, with the datatypes package loaded as the toolkit itself needs. It
reads a standard ISA `metadata.csv`, takes its `Source` and `algo_*`
columns, and copies PRELIM's labelling code — NaN handling, the zero-score
`eps` substitution and all — so what it reports is what PRELIM will
actually produce. On six MATILDA problem datasets its labels match the
toolkit's own `PRELIM` in every cell, at 202 candidate thresholds, in
MATLAB R2026a and in Octave 11.3.

```matlab
% absolute: values are in the metric's own units (here, seconds)
threshold_sweep('runtimes.csv', 'MaxPerf', false, 'AbsPerf', true, ...
    'Epsilon', [1 10 60 600 3600])

% relative: values are fractions of the per-instance best
threshold_sweep('runtimes.csv', 'MaxPerf', false, 'AbsPerf', false, ...
    'Epsilon', [0 0.01 0.02 0.05 0.10])
```

Called without an output it prints the report; `r = threshold_sweep(...)`
returns the same results as a struct instead.
`references/sweep-columns.md` describes every field. What follows is the
part you need whether or not you load it.

**Read `r.data` first.** It is computed before any threshold and is
identical for absolute or relative, which is what lets it rule out a mode
before you commit to one. Three of its fields settle things on their own:

- `ybestZeroFrac` above 0.05 **kills relative mode**; PRELIM raises
  `ISA:PRELIM:manyZeroBest` at the same cut.
- `hasNegative` or `crossesZero` makes the relative ratio meaningless.
- `safeRungInterval` and `pointMassNeighbourBelow` give the range an
  absolute value may occupy near a sentinel without swallowing real data.
  Read them before proposing anything near one; `cutoff - 1` is only safe
  when that interval says so.

**`screenScore` is a screening heuristic, not an objective.** It is the
geometric mean of balance, informativeness and complementarity, and it
exists so that a mechanical search over values has something to maximise.
Use it to find a promising neighbourhood and to throw out degenerate
candidates. Never report it as the reason for a choice, and never tune
anything against it: the justification comes from the argument for the
value, not from this number.

**A high `groupSpread` is not on its own disqualifying.** Benchmark
families are often *designed* to span difficulty, so a label that
separates easy groups from hard ones may simply be correct. The
degeneracy is when the threshold separates groups *instead of*
algorithms. Read the two together:

| `groupSpread` | algorithm discrimination | verdict |
|---|---|---|
| high | strong (wide per-algorithm spread, low `meanJaccard`) | fine — the groups differ in difficulty *and* the algorithms differ |
| high | weak (algorithms failing or succeeding together) | **reject** — the label is a class detector |
| low | either | no group confound to worry about |

A high `|groupRho|` strengthens the rejection, because it means the label
tracks the suite's own design axis; `groupRho` near zero with a high
spread usually means the groups differ but not monotonically, which is
benign. **With more than one family present, check per-family good-rates
too** — a value that is healthy pooled can be degenerate on the smaller
family, and the per-algorithm ranking can differ between families in a
way pooled numbers hide entirely.

**A flat curve at a saturated good-rate is not a plateau.** If 97% of
pairs are already labelled good, the labels have stopped changing and of
course the curve is flat — that is degeneracy, not robustness. The script
returns `plateauWidth` as `NaN` in that range rather than a width that
would be read as evidence.

### Read the two against each other

**Test both values either side, never just the nearer one.** A sweep
optimum at 2834 s sits between 1800 s and 3600 s; evaluate both. This
produces the plateau evidence as a by-product:

- **Similar diagnostics** → plateau. The exact value demonstrably is not
  doing the work. This is the strongest available defence, stronger than
  the justification alone. Report it.
- **Sharply different** → cliff. Neither is safe; see
  `references/failure-modes.md`.

How far to trust the domain proposal tracks whether a hard anchor exists.
With one, it is strong and defensible without the data. Without one it
rests on a community convention, and **that is where it fails** —
fluently and with confidence, which is what makes it hard to catch. A
convention imported from the literature carries the instance population
it was calibrated on, and on a differently-shaped suite it can partition
instances instead of algorithms. So with no hard anchor, weight the
diagnostics much more heavily and treat the proposal as a hypothesis.

**If they agree**, report it plainly: the conventional threshold is also
the one that best discriminates this portfolio. Done.

**If they disagree**, that is data, not a problem to resolve by splitting
the difference. Three readings, to be distinguished rather than averaged:

1. **The portfolio has outgrown the convention.** Everything is now
   within 1%, so 1% no longer separates anything and 0.1% is the live
   threshold. A finding about the state of the art.
2. **The convention is valid but does not transfer** — your instance set
   and the one that established it are not the same population, so a
   value calibrated there partitions instance *classes* here rather than
   algorithms. The tell is a high `groupSpread`/`groupRho`. The fix is
   the mode or a per-instance bound, never a different `epsilon`.
3. **The convention is stale** or was never well-founded. A finding about
   the field. Distinguish this from (2): there the convention is sound
   and merely out of population; here it was never sound.

### When something is wrong

Each of these has its diagnosis and remedy in
`references/failure-modes.md`. Read the entry that matches what you are
seeing.

- Every value is either degenerate or saturated, with nothing in between.
- The scale reaches or crosses zero.
- The metric is censored — a timeout, a cutoff, a penalised score.
- The metric is integer-valued, you are in relative mode, and the
  smallest instances are small enough that a few percent of the best
  rounds to nothing.
- The question is hardness and you have no bound outside the portfolio.

---

## 5. Report it

**Lead with the answer.** The reader wants the setting, not the
derivation. Everything above is working material and almost none of it
belongs at the top — no "I began by considering", and no algebra above
the evidence table.

**Several candidates are welcome**, since the reader is usually the
domain expert and is better placed than you to weigh a trade-off you can
only partly judge. Hedging is not: two values resting on one argument are
a sensitivity range, and belong in the table as one. Three
candidates is the practical ceiling — past that you have stopped choosing
and started surveying.

**1. The candidates.** For each: a config block, one plain
sentence, then at most four bullets — the argument for it, **what the
resulting space will let the author claim and what it rules out**, and
what would change it.

```matlab
opts.perf.MaxPerf = false;
opts.perf.AbsPerf = true;
opts.perf.epsilon = 3599;      % illustrative: a 3600 s budget
```

> An algorithm is good if its runtime is under 3599 s — that is, if it
> proves optimality inside the 3600 s budget.

Phrase that sentence the way PRELIM phrases it in its own log ("An
algorithm is good if its performance is less than X", "... within X% of
the best"), so the reader can check the recommendation against what the
pipeline prints when it runs. Give the encoded value and the domain
meaning together where they differ, as above: the reader needs 3599 to
check the log and 3600 to judge whether the threshold is the right one.

**2. Evidence.** One table covering every candidate and the values either
side of it, to three significant figures, with the candidates marked.
Head the columns in words rather than the script's identifiers —
`overlap`, not `meanJaccard` — and map them to `r.sweep`'s names in one
line underneath. Name every diagnostic in words and say which way is
good: not "meanJaccard 0.393" but "their good-sets overlap by 39%, so
they are largely complementary". If a diagnostic needs its direction
explained twice, drop the number and keep the sentence.

**3. Choosing between them.** Name the **question each candidate
answers** and the fact that would settle it — say, that the strict value
serves the hard-region claim while the looser one serves selection, and
that a run-status flag would decide between them. State a preference if
you have one; if they are genuinely equivalent and the choice is a domain
judgement, say *that*. What you may not do is leave the reader to infer
the difference.

**4. Caveats to carry into the write-up.** Unverified noise floors,
sentinels that mix mechanisms, group rankings that differ from the pooled
ranking, passengers in the portfolio. This is the part people reuse.

**5. Rejected.** One line each, at most three. "Relative mode: inverts on
censored data — all-fail instances label all-good."

Always give the sensitivity range alongside the candidates. It is a
robustness claim, not a menu: if the space only behaves sensibly at one
exact `epsilon`, the range is what shows it.

Two things never to do. **Do not tune `epsilon` against downstream
results** — footprint purity, PYTHIA accuracy and `epsilon` are all free
parameters, and adjusting one to flatter the others produces a beautiful
space that reports nothing falsifiable. **Do not revise it silently.**
Revisiting after seeing downstream results is legitimate, since ISA is
iterative; doing it without saying so is not, because the reader cannot
then tell a threshold that was chosen from one that was fitted. And if a
value came out of the sweep rather than the domain, say so.

### Reviewing an existing threshold

Open with a one-word verdict and a one-sentence reason, then the same
shape:

> **CHANGE.** Relative mode on a censored metric labels all four models
> good on the 35% of instances where none of them solved anything.

Use **KEEP** when the threshold survives the diagnostics, **CHANGE** when
it is degenerate or inverted, and **INVESTIGATE** when the deciding fact
is missing rather than unfavourable. Do not soften a CHANGE into a list
of considerations.
