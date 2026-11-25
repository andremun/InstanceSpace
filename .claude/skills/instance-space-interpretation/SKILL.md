---
name: instance-space-interpretation
description: >
  Interpret an instance space that has already been built, and write up
  what it says about the problem. Use this skill when the user has a
  finished model (model.mat, output CSVs, or the PNG plots) and wants to
  know what it means: what the axes represent, where each algorithm is
  strong or weak, which features explain performance, whether the
  benchmark suite has gaps, and whether an algorithm selector
  outperforms simply running one algorithm everywhere. Trigger on
  "interpret my instance space", "what does this space tell me", "read
  these results", "what do the axes mean", "which features explain
  performance", "is the selector any good", or a request to write up ISA
  results.
  For running the pipeline, choosing options, or the methodology behind
  each stage, use the `instance-space-analysis` skill instead. For
  setting or justifying the good-performance threshold, use
  `performance-thresholds`. For what to change about a space once it has
  been read -- gaps worth filling, a feature set or portfolio that is the
  real limitation -- use `instance-space-improvement`.
---

# Interpreting an Instance Space

The space is already built. The job now is to say what it means, in
language a domain expert can act on.

This skill assumes the `andremun/InstanceSpace` MATLAB toolkit. Re-check
the field names and output filenames below if the toolkit has moved on.

---

## Writing rules — these apply to everything you produce

The reader is a domain expert in their own problem, not in this
toolkit. Write for them.

- **Explain every term the first time.** Not "the beta-hard footprint"
  but "the region where most algorithms in the portfolio fail".
- **Give the consequence, not the setting.** Not "epsilon was 0.05 with
  a relative measure" but "an algorithm counted as good if it came
  within 5% of the best algorithm on that instance".
- **Keep symbols out of the summary.** Axis names, matrix names and
  option fields belong in the detail sections, not in the opening
  paragraphs. Say "the horizontal axis" before you say "z1", and only
  use "z1" afterwards if it saves real space.
- **Quote feature values in their original units** whenever you make a
  claim about a range. A claim expressed in rescaled units cannot be
  acted on.
- **Say when you don't know.** An honest "the features do not explain
  this region" is a finding. A confident story built on evidence that
  is not there is not.

---

## Step 0 — Gather the inputs

Ask once, up front, for everything missing — the required inputs and
the optional context both, since it is the optional context that turns a
description of the geometry into an interpretation. **When something is
missing** below says how to put the request. Wait for that reply, then
proceed with what you have and record what is still missing in the
limitations section. Do not stall beyond that one round of questions.

### Required

| Input | Why it is needed |
|---|---|
| `model.mat` (or the output CSVs and PNGs) | Everything else is read from it |
| Problem description | Without it you can describe geometry but not meaning |
| Performance measure semantics | What was measured, in what units, and **whether runs were cut off by a time or iteration limit** |

The cut-off question matters more than it looks. If runs were capped,
then "bad performance" mixes together *slow* and *never finished*, and a
boundary in the picture may be tracing the cap rather than anything
about the algorithms.

### Improves the interpretation a lot

- **How the instances were obtained or generated.** Which generators,
  which parameter ranges, which real-world collection.
- **The console log from the build run.** Several things are printed and
  never saved: the percentage of instances where several algorithms tied
  for best and the tie was broken at random, features dropped for having
  too many missing values, a warning that the feature clustering was
  poor, and a warning that the outer boundary fell back to a simpler
  calculation. Each changes how much weight a result deserves.
- **What the algorithms are, and how they relate.** Two variants of one
  solver sharing a region is expected. Two unrelated methods sharing one
  is a finding. Also ask whether they were tuned or run at defaults — an
  untuned algorithm's region describes its default settings.
- **Why the good-performance threshold was set where it was.** The value
  is in the model but in none of the CSVs; the reasoning is in neither,
  and every region in the picture depends on it.
- **Feature meanings** for any feature whose name does not make its
  meaning obvious — or the feature generation code, which settles a
  definition instead of describing it. With neither, and the raw instance
  files to hand, the `feature-recovery` skill recovers a definition from
  the values.
- **The toolkit version or commit.** The model records no version field.
  Feature selection, the classifier stage and the footprint stage have
  all changed behaviour across recent versions, and so have the success
  rates in the summary table (Steps 2 and 8).

### Worth asking, cheap to answer

- **What the interpretation is for** — a paper's results section, a
  referee reply, a decision about adding an algorithm, or a defence of
  the benchmark suite. Same space, different emphasis.
- **What the user expected to see**, and any published finding this
  should agree or disagree with. A contradiction is then reported as a
  finding rather than smoothed over — but first establish that the model
  *is* the run the write-up describes, because a stored model and its
  published account can diverge. Four cheap checks settle it: the
  feature count in the metadata header against the published feature
  table, the instance count, the classifier (Step 8), and
  `model.pilot.A` against any published projection matrix. If any of
  them disagree, later differences are version differences and must not
  be reported as contradictions.
- **Any `explore()` results** on a separate test set. Much stronger
  evidence than the training instances alone.

### Already there — never ask for these

**With `model.mat`:** all settings (`model.opts`), instance names,
feature names, algorithm names, the raw and processed feature values,
raw performance, the good/bad labels, the best algorithm per instance,
and every stage's output. The CSVs and PNGs are convenience only;
`obj.plot(...)` regenerates the figures from the model.

**With only the exported CSVs and PNGs**, nearly all of that survives —
*Where everything lives* names the file for each one. Two things the user
would have to supply do not, so both move onto the list you ask about:

- **The settings.** `model.opts` is written to no CSV, so the
  good-performance threshold's kind and value have to be asked for, not
  only the reasoning behind them — and every region in the picture rests
  on that threshold. `options.json` settles it if the run kept it.
- **Which classifier was trained** — see *The `svm` in those names is a
  default, not a guarantee* under *Where everything lives*.

`obj.plot(...)` is unavailable too, so you are limited to whichever PNGs
were written.

### When something is missing

**Tell the user what is missing and ask for it — once — before you
produce the write-up.** A missing input you notice while reading the
model is the trigger for that ask, not something to route straight into
the limitations section.

Put every missing input you have noticed into a single message: name what
you do not have, say in a few words what each one would let you settle,
and say that you will proceed either way. Then stop, and let the user
answer. This is the one place in the skill where it is right to wait, and
one round is all it takes — whatever comes back, the next reply starts
the write-up.

Fold in whatever arrives. *"I don't have that"* and *"just get on with
it"* are complete answers, not failures: they cost the user a line, they
tell you the input will not arrive rather than having been overlooked,
and they are the expected reply to at least half the list. So do not ask
twice, do not re-ask the same thing in a different form later, and never
make a second round a condition of delivering the interpretation.

The failure this guards against is a quiet one. A folder holding a
model, a problem description and every output CSV looks complete, so
nothing ever feels missing enough to be worth raising — and the write-up
ships with a limitations section listing context the user could have
supplied in a sentence, had anyone asked.

One case must be called out explicitly rather than worked around: if the
metadata had no `Source` column, `model.data.S` does not exist, and
**Step 3's source analysis cannot be done at all**. Say so plainly
instead of substituting something weaker.

---

## Where everything lives

| What you want | In the model | In an export |
|---|---|---|
| Settings actually used | `model.opts` | — (`options.json`, if the run kept it) |
| Instance / feature / algorithm names | `model.data.instlabels` / `.featlabels` / `.algolabels` | the row and column headers of any export |
| Raw feature values | `model.data.Xraw` | `feature_raw.csv` |
| Feature values used by the pipeline | `model.data.X` (trimmed and rescaled, not raw) | `feature_process.csv` |
| Raw performance | `model.data.Yraw` | `algorithm_raw.csv` |
| Good/bad labels | `model.data.Ybin` | `algorithm_bin.csv` |
| Best algorithm per instance | `model.data.P` | `portfolio.csv` |
| How many algorithms are good per instance | `model.data.numGoodAlgos` | `good_algos.csv` |
| Instance source/class | `model.data.S` (only if the CSV had a `Source` column) | the metadata's `Source` column |
| Which features survived selection | `model.sifted.selvars` | — (infer from `feature_process.csv` headers) |
| Suggested better cluster count | `model.sifted.Ksuggested` (written only when the clustering that ran scored poorly — also absent when the stage was skipped) | — |
| Axis definitions | `model.pilot.A` | `projection_matrix.csv` |
| Instance coordinates | `model.pilot.Z` | `coordinates.csv` |
| Performance loadings | `model.pilot.C` | — |
| Outer boundary | `model.cloist.Zecorr` (correlation-pruned, prefer it) or `.Zedge` (feature ranges alone, wider) | `bounds_prunned.csv` / `bounds.csv` |
| Per-algorithm and Selector/Oracle table | `model.pythia.summary` | `classifier_table.csv` |
| Selector accuracy, precision, recall | `model.pythia.accuracy` / `.precision` / `.recall` | the `CV_model_*` columns of `classifier_table.csv` |
| Per-algorithm predicted-good labels | `model.pythia.Yhat` | `algorithm_svm.csv` |
| Selector's pick per instance | `model.pythia.selection0` (blank when unsure), `.selection1` (always picks) | `portfolio_svm.csv` gives `selection0`; rebuild `selection1` by filling its blanks with the highest-mean column of `algorithm_bin.csv` |
| Footprint measurements | `model.trace.summary` | `footprint_performance.csv` |
| Footprints per algorithm | `model.trace.good{i}`, `.best{i}` | `footprint_<algo>_good.csv` / `_best.csv` |
| Region where few algorithms are good | `model.trace.hard` — at most `opts.perf.betaThreshold` of the portfolio, 55% by default, so not the all-fail region | — (`beta_easy.csv` gives the easy/hard split per instance) |
| Every projection attempt, and its score | `model.pilot.alpha` / `.perf` | — |
| New instances projected into this space | `obj.getResults(idx)` — projecting does not alter `model` | — |

The dashes matter: six of these exist only in the model, so an export
cannot answer what the settings were, which features selection dropped,
or how the projection scored against its own rejected attempts.

Figures written by the pipeline, if `outputs.png` was on:

| File | Shows |
|---|---|
| `distribution_feature_<name>.png` | One feature's values across the picture |
| `distribution_performance_global_normalized_<algo>.png` | Performance on one shared scale across algorithms |
| `distribution_performance_individual_normalized_<algo>.png` | Performance rescaled per algorithm |
| `binary_performance_<algo>.png` | Where that algorithm is good |
| `binary_classifier_<algo>.png` | Where the selector *predicts* it is good |
| `footprint_<algo>.png` | Its footprint against its actual good instances |
| `distribution_number_good_algos.png` | How many algorithms succeed per instance — the hard/easy map |
| `distribution_beta_score.png` | Hard/easy as a yes/no split |
| `distribution_portfolio.png` | Which algorithm is actually best per instance |
| `distribution_svm_portfolio.png` | Which algorithm the selector picks — drawn from `selection0`, so blank where it abstained |
| `distribution_boundary.png` | The outer boundary (2D projections only) |
| `distribution_sources.png` | Instance sources |

Interactively: `obj.plot('sources')`, `'portfolio'`, `'good'`,
`'footprint'`, `'boundary'`; the last two take an algorithm index.

**The `svm` in those names is a default, not a guarantee.** PYTHIA writes
`distribution_svm_portfolio.png`, `algorithm_svm.csv` and
`portfolio_svm.csv` under those names whatever classifier it trained.
Read them at face value — assume an SVM — unless
something that records the run says otherwise: `model.pythia.summary`
(Step 8), `options.json`, or the build log. Those override the
filenames. With only the exported CSVs and PNGs there is nothing to
override them, so treat the SVM as an assumption.

**Do not build claims on the best-algorithm outputs.** PRELIM breaks a
tie for best algorithm by picking one of the tied algorithms at random,
so `distribution_portfolio.png`, `footprint_portfolio.png`,
`portfolio.csv` and every `_best` footprint depend on an RNG draw and can
change between runs on the same data. Nothing good/bad-based is affected,
because a tie means every tied algorithm counts as good: use
`binary_performance_<algo>.png`, `algorithm_bin.csv`,
`distribution_number_good_algos.png` and the `_good` footprints instead.
Where a question genuinely needs "which algorithm is best here", answer
it from the good/bad maps and the performance values, and say ties were
broken arbitrarily.

**Two colour-scale traps.** The shared-scale performance plots apply a
logarithmic transform before colouring, which suits runtimes but
compresses and can mislead for a measure that is already a bounded
score. The per-algorithm plots rescale each algorithm separately, so
colours **cannot be compared between algorithms** — a dark region in one
plot and a dark region in another do not mean the same performance. And
the feature plots colour the trimmed, rescaled feature values, not the
original units.

### Pipeline defaults

Any option a run did not set takes its default, whether it was left
out of `options.json` deliberately or simply never thought about. If
you have the model, you do not need this table for that run:
`model.opts` records the value every option actually took, defaults
included. Without the model, `options.json` tells you only the options
it names, and everything it leaves out took its default. That is the
case this table is for.

One catch on the model side: when the toolkit loads an older model, it
fills in any option the model does not have with *today's* default. On
a model saved before an option existed, `model.opts` can show a value
the run never used.

**This table is a snapshot, checked against the toolkit on
2026-11-23.** The source of truth is `utils/ISAdefaults.m`, one line
per option. If a conclusion depends on a default, check the value
there. The file can be older or newer than the version the user ran,
so if you still cannot see what their run used, ask them.

| Option | Default |
|---|---|
| `perf.MaxPerf` | false (minimise) |
| `perf.AbsPerf` | false (relative) |
| `perf.epsilon` | 0.05 |
| `perf.betaThreshold` | 0.55 |
| `prelim.nanThreshold` | 0.20 |
| `prelim.iqrMultiplier` | 5 |
| `sifted.rho` | 0.10 |
| `sifted.pval` | 0.05 |
| `sifted.K` | 10 |
| `sifted.diagnostics` | true |
| `cloister.maxFeatures` | 20 |
| `trace.method` | `trace3` |
| `trace.PI` (purity) | 0.6 |
| `trace.minInstances` | 4 |
| `trace.minAreaFrac` | 0.01 |

**The skills also quote some of these values where an argument depends
on the number**, usually as "N by default". If the file changes,
update this table and its date, then search the skills for "by
default" and fix each copy. Watch for numbers written as words, such
as "come out as ten".

---

## Step 1 — Describe what was analysed

Be detailed. This section is what makes the rest readable.

- **Instances**: how many, where they came from, what a single instance
  is in this problem.
- **Features**: how many were supplied, and what each one measures in
  plain words. Where a feature's meaning is unknown and cannot be
  inferred from its name, say so rather than guessing.
- **Algorithms**: how many, what each is, and which ones are related to
  each other.
- **Performance measure**: what was recorded, in what units, whether
  higher or lower is better, and whether runs were cut off.
- **The definition of good performance, in a sentence.** Either "good
  means beating a fixed threshold of X" or "good means coming within X%
  of the best algorithm on that same instance". Then note the
  consequence: with a relative definition, almost every instance has at
  least one good algorithm by construction, so "hard" later on means
  "few algorithms do well here", not "nothing works".

Close with a short note on what actually ran, because low and high
feature counts both change the pipeline quietly:

- The feature-selection stage skips its search and keeps whatever the
  correlation filter left whenever the survivors do not outnumber
  `opts.sifted.K`, 10 by default. That is the gate that usually fires;
  three or fewer features skips it too, but rarely applies.
- Above the feature cap (`opts.cloister.maxFeatures`, 20 by default)
  CLOISTER cannot enumerate the range corners, so the outer boundary becomes the
  plain convex hull of the projected instances and a
  `tooManyFeatures` warning fires. Below the cap `Zedge` comes from the
  observed feature ranges and `Zecorr` keeps only the corners consistent
  with the pairwise correlations — both wider than the instances,
  neither a theoretical limit. Every "the space is under-filled" claim
  depends on knowing which boundary you are measuring against.
- The random tie-breaking rate is measured in two places and the two
  legitimately differ: the log reports it before the density filter runs,
  while the model holds only what survived the filter. If the filter
  removed most of the tied instances, say which it discarded — that is a
  finding, not a caveat.

---

## Step 2 — How the algorithms do overall

Read `model.pythia.summary`. It has one row per algorithm plus an
**Oracle** row (perfect choices) and a **Selector** row. Report, per
algorithm: the proportion of instances it is good on
(`Probability_of_good`), and its average performance and spread.

That proportion summarises `model.data.Ybin` (`algorithm_bin.csv`), a
good/bad label for every instance against every algorithm. Count from
that matrix for exact numerators and denominators, or for any breakdown
the summary table does not carry. It cannot separate a missing run from a
bad one: PRELIM substitutes ±Inf for a missing value before thresholding,
so both land as a zero. The raw performance (`model.data.Y`,
`algorithm_raw.csv`) is where absent entries show.

If any performance values are missing, check the denominator of that
proportion. Before v0.9.2 the toolkit counted a missing value as *not
good*, which lowers the rate of every algorithm that has any. From
v0.9.2 it divides by the observed instances only. Recompute over the
observed instances whenever the performance data is incomplete.

This is the baseline for everything that follows. An algorithm good on
most instances is *expected* to cover much of the picture; a region of
the same size means far more for one that is rarely good.

Name the best and worst performers overall here. Step 6 then asks
whether those overall rankings actually hold everywhere.

---

## Step 3 — Where the instances sit

Look at `distribution_sources.png`, or `obj.plot('sources')`.

- Do real-world instances land somewhere different from generated ones?
- Do randomly generated instances bunch near the middle? This is the
  usual pattern and it means the generator produces "average" instances
  rather than diverse ones.
- Does any one source occupy a distinctive region? That region's
  features are what makes that source distinctive — connect it to
  Step 5.
- Are there sizeable empty areas inside the outer boundary? List them by
  location. They go in the write-up (Step 9), and they are where
  `instance-space-improvement` starts.

If there is no source information, skip this and say so.

---

## Step 4 — Which features survived, and what the axes mean

**Feature selection.** If the pipeline reduced the feature set, list
what survived and what did not. A dropped feature is not necessarily
unimportant; Step 5 says what to do about that. If `model.sifted.Ksuggested`
exists, the feature grouping was poor and a different grouping was
recommended; mention it, because it weakens everything downstream.

To list *dropped* features you need the original feature names, from the
metadata CSV headers — the model's `featlabels` holds the survivors only.

**The axes.** `model.pilot.A` gives each axis a weight for every
surviving feature. Each axis is a weighted blend of features. For each
axis, rank the features by the size of their weight and translate the
top few into a sentence:

> Instances at the right of the space are higher on feature B and lower
> on feature C.

That sentence is the foundation of the whole report. Everything after it
refers back.

Three cautions, all of which belong in the write-up:

- The weights apply to features **after** outlier trimming and
  rescaling, so they describe the adjusted values, not raw units.
- **The axis directions are not fixed by the method.** The fit is
  seeded, so repeating it on the same data reproduces the same picture
  — but changing the instance set can flip left for right, and so can a
  different number of fitting attempts. Always confirm the direction
  against an actual feature plot before writing a directional claim —
  never read it from the sign of a weight alone.
- The projection is a compromise fit. It cannot show every feature
  clearly, which is exactly what Step 5 looks for.

Also worth reading here: `model.pilot.C` links position in the picture
to algorithm performance, which gives each algorithm a rough direction
of improvement. Treat it as a cross-check on Step 6 rather than as the
primary evidence, and note that performance is converted to a shortfall
from the best algorithm before this is fitted, so its directions point
the opposite way to raw performance.

---

## Step 5 — How each feature spreads across the picture

**Triage first when there are many features.** Ten or twenty feature
plots is a lot to eyeball, and the muddy ones are hard to spot by
comparison when you view them one at a time. Correlate each feature's
values against the two coordinates first; that ranks the features and
tells you which plots deserve a careful look and which are candidates
for the muddy pile. Then confirm visually — the correlation is triage,
not the finding, because a strong non-linear pattern (highest in the
centre, lowest at both edges) scores near zero on a correlation and is
still a perfectly clear gradient. Keep the numbers out of the write-up;
the report describes what the plots show.

Go through `distribution_feature_<name>.png` one at a time. For each
feature, describe its pattern in plain directional language:

- "increases steadily from left to right"
- "decreases from the top-right corner towards the bottom-left"
- "highest in the centre, lower at the edges"
- "no clear pattern — high and low values are mixed throughout"

Then sort the features into two groups:

- **Clear features** — a visible, consistent gradient across the
  picture. These are the ones you can use to describe regions, and the
  ones that make claims in Step 7 meaningful.
- **Muddy features** — no consistent pattern, high and low values
  mixed throughout.

The projection stage tries to make features vary smoothly across the
picture, so a muddy feature is one it failed on. That is worth
reporting, not hiding: it usually means the feature is unrelated to the
others, or measures something two axes cannot show at once. If most
features are muddy, the projection fits the data poorly and the whole
interpretation should be presented with that caveat.

**Correlate every muddy feature against performance before setting it
aside, and every feature selection dropped along with them.** "Muddy"
describes the projection, not the feature: it means two axes could not
show this feature clearly, which is not the same as the feature not
mattering. One correlation per (feature, algorithm) pair against the
binary labels settles it in seconds. Where one comes back strong, report
it as a finding, flagged as resting on the raw values rather than on the
picture — this is one of the few places the map actively misleads, since
a feature with no gradient to point at looks like a feature that does
not matter.

---

## Drawing your own plots

The shipped PNGs may not exist. Draw your own from `coordinates.csv` (or
`model.pilot.Z`) plus whatever you are colouring by — the projected
coordinates are just two columns, and every map in this skill is a
scatter of those two columns with a different colour rule.

Do it even when the PNGs do exist, if a different view would answer the
question better. The shipped set is one plot per feature and one per
algorithm; it has no map of which algorithm is *uniquely* able to solve
an instance, no map of what the selector recommends against what it
should recommend, and no view that marks the regions you identified in
Step 6. Those are often the most informative pictures in the whole
report.

Keep them honest: same projection and same axis orientation across every
panel, so they can be read against each other.

**Match the figure to the performance measure.** A scatter of times is
the wrong figure for a censored runtime (everything that did not finish
recorded at the cutoff), because the cutoff pile dominates. The right
figure is the one that field already uses: cumulative instances solved
against the time budget, one line per algorithm, log time axis — a
cactus plot. It shows at a glance whether a curve has flattened, which
tells you whether more time would have solved more instances, and it
exposes things a good/bad table cannot.

---

## Step 6 — How performance spreads across the picture

Use `binary_performance_<algo>.png` for good/bad, the shared-scale
performance plots to compare algorithms, and
`distribution_number_good_algos.png` for difficulty.

**Per algorithm**: where is it strong, where weak, stated as regions of
the picture ("strong across the whole left half", "weak in the top-right
corner"). Then read those against Step 2's ranking, which is an average
and hides exactly this — where do strong algorithms fail and what covers
them, where do weak ones win alone. An algorithm with no region of its
own is a finding too: a candidate for dropping from the portfolio, and a
reason the selector may have nothing to gain (Step 8).

**Difficulty**: `distribution_number_good_algos.png` shows how many
algorithms succeed on each instance, so dark regions are where the whole
portfolio struggles. Is there a clear hard region and a clear easy
region? Where do they sit relative to the feature gradients from Step 5?
This is often the single most interesting finding in a space.

Note the caveat: this reads as genuine difficulty only when good
performance was defined against a fixed threshold. With a relative
definition, every instance has a good algorithm by construction, so the
map shows how *many* algorithms cope, not whether the instance is hard
in absolute terms.

**Footprints.** Read `model.trace.summary` alongside the plots
(`footprint_performance.csv` carries these columns in an export). Each
algorithm gets two regions there, a *good* footprint and a *best* one —
use the good one, for the reason given under *Where everything lives*.
It comes with three measurements:

- **Area** — how much of the picture it covers. Breadth of competence.
- **Density** — how tightly instances are packed in it, relative to the
  picture overall. Above 1 means denser than average, so there is real
  evidence behind it. At or below 1 the region is thinly supported.
- **Purity** — the share of instances inside it the algorithm is
  actually good on. Below roughly 0.7 the region is contradicted often
  enough that claims about it should be softened.

A region that is large but thin and impure is weak evidence, not a
strong result. And read every area against that algorithm's overall
success rate from Step 2, never on its own.

**Overlaps** between algorithms are informative, not a failure. Two
related solvers overlapping is expected. Two unrelated methods
overlapping heavily means they are genuinely competitive across the same
instances, which is itself the finding.

**A missing footprint is not proof of a weak algorithm.** It can be an
artefact of the footprint calculation, and it is the expected result for
an algorithm that was not in the trained model but appears in new data.
Check before concluding.

---

## Step 7 — Link features to performance

This is where the "why" comes from. Combine Steps 5 and 6 into claims
specific enough that they could be wrong:

> Algorithm A is good when feature B is above *x* and feature C is
> below *y*.

Rules for these claims:

- Build claims *about the picture* only on **clear** features from
  Step 5. A muddy feature, or one that selection dropped, cannot support
  a claim about the picture — but consult Step 5's correlations before
  dropping the claim, since the raw values may support it well. Where
  they do, keep it and say which evidence carries it. Of the two ways to
  be wrong here — keeping a claim the picture does not show, or dropping
  an effect that is really there — the second is worse.
- Quote ranges in **original units**. Recover them from the raw values
  for the instances in that region, not from the rescaled values used
  internally.
- Where no clean feature story exists, fall back to a location claim:
  "algorithm A is good in the lower half of the picture, below about
  -2 on the vertical axis". Weaker, but honest, and still usable for
  selection.
- Where neither works, say the features do not explain that region. That
  points at a missing feature and is a real result.

**Checking a claim before you commit to it (recommended, not required).**
If you have the model to hand these are quick, and each has caught a
false claim before:

- Pull out the instances the claim covers and confirm the performance
  gap is actually there in the raw numbers.
- Check the claim is not carried by a single instance source or one
  tight cluster of near-duplicates.
- Check it survives if the region is drawn slightly differently.
- Check it is not an artefact of ties being broken at random, or of a
  time cut-off.

A claim that fails a check should be reported as unsupported rather than
quietly dropped.

---

## Step 8 — How does the selector perform?

### Start with the per-algorithm models

Each algorithm has a classification model that predicts whether that
algorithm will be good on an instance, from where the instance sits in
the space.

Present the classifier results for each algorithm in one figure — a panel
per algorithm from `binary_classifier_<algo>.png` — with a few brief
comments on where each is expected to do well. Then give the statistics
behind them: accuracy, precision and recall per algorithm
(`model.pythia.accuracy`, `.precision`, `.recall`). These three are
cross-validated, unlike everything else the step reports, so they are the
closest thing here to an out-of-sample estimate.

**Which classifier these are.** Selectable, and defaulted to
k-nearest-neighbour in recent versions, so "the selector" is not always an
SVM. The last two columns of `model.pythia.summary` name it:
`NumNeighbors` / `Distance` for k-nearest-neighbour, `BoxConstraint` /
`KernelScale` for an SVM. If the summary is unavailable, the build log
names it outright.

### How they become one selector

If at least one algorithm is predicted good on an instance, the selector
picks whichever of those has the highest cross-validated precision. So
precision is not only a score — it decides who wins whenever more than
one model says good, and a weak model demotes its algorithm everywhere
rather than only where it is wrong. Where no model predicts good there is
nothing to pick, which is why there is an abstention rate and two
outputs:

- `selection0` — the *confident* output. Blank wherever no model
  predicted good.
- `selection1` — the *operational* output. Always names something,
  falling back to the algorithm most often good across the training set.

The summary table's `Probability_of_good` for the Selector row uses
`selection1`, so it is an always-picks figure and compares fairly against
always running one algorithm, or against any rule that always chooses.
Its denominator is every instance, the fallback's hits and misses
included. So do not put it beside the abstention rate without saying
which is which: side by side, the two read as though the score had been
earned only on the instances the selector was willing to answer, which
would make it a much better selector than it is.

### What to compare it against

- **each algorithm run everywhere**, and in particular the **best single
  algorithm** — that one is the bar,
- the **Oracle** — perfect choices every time, the ceiling.

Score them all the same way, and score them twice: on the good/bad rate
the toolkit reports, and on the performance measure itself. Those two can
disagree, and where they do the second is usually the one the user cares
about.

**Recompute the Oracle's success rate before using it.** Versions before
v0.9.2 set the Oracle's `Probability_of_good` to 1 in the summary table
instead of working it out from the data. v0.9.2 computes it, but the model
records no version field, so a stored model cannot say which kind it is.
Recompute it every time; on a current model the two agree. The literal 1
is right only when good performance was defined *relative* to the best
algorithm on each instance, where the best algorithm always clears the bar
by construction. Under an **absolute** bar it is wrong: instances where
even the best algorithm misses the threshold have no good algorithm at
all, and the Oracle cannot invent one. Compute it as the fraction of
instances with at least one good algorithm (`mean(any(Ybin,2))`) and use
that throughout.

This is not cosmetic. It inflates the ceiling, so the headroom above the
best single algorithm looks far larger than it is and the selector looks
worse than it is. A quick internal check catches it: the Oracle's average
performance is worked out correctly from the per-instance best, so if
that average sits worse than the good-performance bar for some instances,
a 100% success rate is impossible on its face.

### Then answer these

1. **Does it beat the best single algorithm?** Give the margin in
   instances as well as percentage points: two rates a fraction of a point
   apart say far less than the count of instances that fraction stands
   for. If it barely beats it, say so — running that one algorithm
   everywhere may be the better recommendation.
2. **How far does it get towards the Oracle?** The best single algorithm
   is the starting point and the Oracle is the ceiling, so where the
   Selector sits between them says how much of the possible improvement
   over running one algorithm everywhere it actually achieves.
3. **Which instances does it get wrong?** Find them on the map and
   describe them as Step 6 described algorithm regions. Separate two
   kinds: instances it abstained on, where the fallback took over
   (`selection0` is blank), and instances it committed to and still got
   wrong.

If it abstains a lot, check the stored hyperparameters
(`model.pythia.param1`, `.param2`) before concluding the features are at
fault. Two things point at the tuning instead: the same tuned values for
algorithms that behave very differently, and a value sitting at the edge
of its search range.

### Last, try writing the selector as a rule

Only once the comparisons above are done. Read the selector's own picks
from `model.pythia.selection1` — the operational output, the one a rule
would be replacing — and ask whether a rule reproduces them, with Step
7's feature and location claims as the obvious candidates: thresholds on
named features, or boundaries drawn on the picture. An export gives
`selection0` instead, in `portfolio_svm.csv`, but `selection1` rebuilds
from it exactly: fill every blank with the algorithm good on the largest
share of instances, which is the highest-mean column of
`algorithm_bin.csv`. That is what the toolkit itself does.

It does not have to fit in one line, and with several algorithms and
several features it rarely will. The test is whether someone could apply
it by hand: a condition or two per algorithm is still a rule, while a tree
with a leaf for every region is a second model and gains you nothing.
Report how often the rule and the selector agree, and score the rule the
same way you scored everything else. A rule fitted to these instances and
scored on them is an in-sample fit — but so is the selector's own success
rate, which comes from the final models predicting the data they were
trained on, so the two are comparable as they stand. Say that, rather
than cross-validating one side of it. An out-of-sample answer needs
held-out instances.

A selector a simple rule reproduces *is* that rule, and the rule is the
more useful thing to hand over: someone can state it, check it and act on
it without ever loading the model. Where nothing simple comes close, say
so — the selector is using the geometry in a way no stateable rule
captures, which is a result in its own right.

---

## Step 9 — Write it up

Suggested structure:

1. **The problem** (Step 0) — what it is, and what one instance
   consists of.
2. **What was analysed** (Step 1) — a short passage per component of the
   metadata:
   - **Instances** — how many, and where they came from.
   - **Features** — how many, and what each measures in plain words.
   - **Algorithms** — how many, what each is, which are related to each
     other.
   - **Performance** — what was measured, in what units, which
     direction is better, and whether runs were cut off.
   - **Good performance** — the definition in one sentence, with its
     consequence.
3. **Headline findings** — two to four, usually as bullet points, each
   specific enough that it could be wrong. No symbols, no toolkit
   vocabulary.
4. **Where the instances sit** (Step 3). Where the sources are known and
   the pattern is worth showing, say where each one falls.
5. **Which features survived, and what the axes mean** (Step 4),
   including the projection matrix itself.
6. **How each feature spreads across the picture** (Step 5).
7. **How performance spreads, and what explains it** (Steps 6 and 7).
8. **How the selector performs** (Step 8).
9. **Limitations** — always present, never optional.

If the user asked what to change rather than only what the space says,
that is `instance-space-improvement`, and it belongs after the
limitations rather than inside them.

The limitations section must cover:

- How much rests on the good-performance definition — what would change
  under a different one.
- Any region the features do not explain.
- Any input that was unavailable, and what could not be concluded
  because of it.

---

## Licence

The `andremun/InstanceSpace` toolkit is distributed under the PolyForm
Noncommercial License 1.0.0 (non-commercial research and educational use
only). This skill documents how to read that toolkit's output and
inherits the same non-commercial spirit.
