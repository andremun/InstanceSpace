---
name: feature-proposal
description: >
  Propose new instance features for Instance Space Analysis, and assess an
  existing feature set for redundancy and gaps. Works cold-start from a
  problem description and instance files, seeded with feature definitions the
  user already has, or — when the job is improving an analysis rather than
  building one from scratch — from a built instance space and its
  interpretation, which say which features were dropped, which ones the
  projection could not show, and which regions nothing explains. Use when
  asked to suggest, propose, invent, or brainstorm features; to review, audit,
  or find gaps in a feature set; to work out what feature is missing for an
  algorithm or a region the current features fail to explain; or to check
  whether a proposed feature is worth computing. Triggers on "propose
  features", "suggest features", "what features should I use", "feature
  engineering", "audit my features", "are my features any good", "what feature
  is missing", "feature gaps", "what feature would explain this region".
---

# Feature proposal

Propose instance features, then measure which hold up. The output is a
feature *specification* — definition, code, reasoning, measured effect — never
a column of numbers on its own.

Three adjacent skills: `feature-recovery` works out what an *existing* column
measures (values decide, names mislead); `performance-thresholds` settles the
good/bad rule this skill needs as an input; and if a built space is what
raised the question, `instance-space-improvement` works out whether the
features are the problem before you get here. For the surrounding methodology
see `instance-space-analysis`.

## Inputs

Ask for anything you are missing before you start — the required inputs and
the optional ones. Send one message, wait for the reply, then get on with it
using whatever you have. *When something is missing*, at the end of this
section, says what that message should look like.

**Required**

- **The problem description** — objective, constraints, what makes an instance
  an instance.
- **The instance files themselves**, plus the format spec.
- **The instance list**, with per-algorithm performance.
- **The good-performance rule** — absolute threshold or relative to best,
  maximise or minimise. If nobody says, propose one and say you have; see
  `performance-thresholds`.

**Optional**

- **What each algorithm does** — a description of its approach, or its source
  code. Usually left out, and the one most worth asking for: without it you
  can only propose generic statistics about an instance's structure, or guess
  at what drives performance from the performance data itself. With it you can
  propose features that measure the structure a given method depends on, and
  those are the ones that make footprints separate rather than just tracking
  difficulty.
- **The existing feature set.** Whether you have this picks the mode below.
- **A built instance space and its interpretation.** For improving an analysis
  rather than building one from scratch. It tells you more than the feature
  definitions do: which features selection dropped, which ones the projection
  could not show, and any region whose performance no feature explains.

### When something is missing

Ask for all of it at once. List what you do not have, say in a few words what
each one would let you do, and make clear you will carry on without it. Then
stop and wait. One round is enough: whatever the answer is, start the work
after it. *"I don't have that"* and *"just get on with it"* are both fine
answers, so do not ask again, and do not slip the same question back in later.

Two things you cannot work around. With no instance files there is nothing to
compute a feature on, and with no per-algorithm performance there is nothing
to correlate against, so Reject and Rank are both out. Writing definitions and
code is still worth doing — just be clear that is what you are handing over,
and that none of it has been measured.

## Two modes

### Seeded (existing features supplied) — the default

Much more efficient. Do this whenever a feature set exists.

1. **Audit the existing set first.** Degenerate and near-constant columns;
   missing-value rates against the 20% cut INIT applies; redundancy clusters, grouped
   the way SIFTED will group them — k-means on correlation distance into
   `opts.sifted.K` clusters, 10 by default, and read from an existing space's
   own `model.opts.sifted` where there is one — and which algorithms each column
   actually explains. Report the clusters, because SIFTED keeps exactly one
   feature from each: a set of thirty columns is going to come out as `K`
   features, ten at the default, whatever you add to it, and which ones
   survive is decided by who shares a cluster with whom. To run this, pass the set itself to `scripts/evaluate.py` as
   `--features` and leave `--existing` off — the degeneracy, correlation and
   clustering reports then describe that set rather than a batch of candidates.
2. **The audit tells you what to look for.** What it finds missing is a far
   better starting point than a hunch: *"propose features that explain
   algorithm C, which nothing here currently does."*
3. **Propose candidates**, then work through the stages below. The overlap
   measurement under Rank applies only in this mode.

**With a built space, the audit comes partly free.** The interpretation
already names the features selection dropped and the ones the projection could
not show, and its Step 7 may have found a region no feature explains at all.
Start from there rather than working it out again: *"propose features that
explain the region nothing currently does"* is sharper than anything the
correlations alone will tell you, because it says where the feature has to
work.

Two cautions when the space is what pointed you here. A feature can correlate with
performance and still not be drawable in the projection;
`instance-space-interpretation` says keep it anyway, but it will not improve
the picture — so say whether you are proposing it to explain performance or to
fix the map. And the space only judged the features it was given: a column
dropped for being redundant may have been redundant with something you are
about to replace.

### Blind (no existing features)

Use this when there is no feature set yet and you have to develop one
yourself.

It is also worth running **once** where features already exist — not to find
new ones, but as a check: if an independent pass arrives at the measures you
already have, both sets look sound. Expect it to confirm rather than discover,
though, so it is not the routine way to find features.

When running it that way, **do not open the existing feature file.** Reading it
throws away the convergence evidence, which is the only thing the pass was for.

## Propose candidates

Breadth first. The checks below are cheap and most candidates will not survive
them, so a short list of carefully polished candidates is the wrong shape to
arrive with. Propose as many as the problem suggests, write them cheaply, and
let the measurements choose between them.

Breadth means different *sources* of structure rather than variations on one
idea. Work down this list and ask what each would be for the problem in front
of you:

- **Scale** — how big the instance is, in whatever units the problem counts in.
- **Distribution shape** — how the quantities in the instance are spread out:
  spiky, flat, clustered, heavy-tailed.
- **Relational structure** — how the elements of the instance constrain one
  another, measured on whatever graph or incidence structure the problem has.
- **Tightness** — how close the instance sits to its binding constraints,
  which is often where easy turns into hard.
- **Heterogeneity** — whether the instance is built from interchangeable parts
  or from a few distinctive ones.
- **Probe results** — run something cheap and measure what it did: the quality
  it reached, how early it stopped improving, how far apart repeated runs
  landed. These cost more to compute than a static statistic, and they measure
  the instance's response to being solved rather than its description.

Where the algorithm descriptions were supplied, aim better than that list
allows: for each algorithm, ask what structure it exploits and what would
break it. Features derived that way are the ones that separate the algorithms
instead of only tracking difficulty, which is what the space needs most and
what a generic statistic rarely delivers.

Stop when new candidates have become restatements of ones you already have,
rather than at a count.

## Reject

Two checks, in this order, both cheap. Neither asks whether a feature is any
good — that is Rank's job. These only remove candidates that cannot work at
all.

1. **Degeneracy.** Drop a constant column — it carries no information — and
   note that at or above `opts.prelim.nanThreshold`, 20% missing by default,
   the pipeline drops the column for you whatever you think of it. A
   near-constant column is a warning, not a rejection: a feature that is 1 on
   exactly the instances one algorithm handles and 0 elsewhere is near-constant
   and perfect. Usually what the imbalance costs is sample size, since the
   correlation rests on the minority rows alone, so take it to the correlation
   check rather than discarding it here.

   One kind of imbalance costs more than that, and a binary feature walks
   straight into it, so check the interquartile range as well. If the middle
   half of the column is a single value then the IQR is zero, median ± 5·IQR
   has zero width, and PRELIM's bounding step clips *every* row to that value
   (`core/PRELIM.m:225`). SIFTED is then looking at a constant column, scores
   it zero (`core/SIFTED.m:117`) and drops it, so the feature never reaches a
   correlation at all. Any binary column whose minority is under a quarter of
   the rows lands here, which includes the one described above. The idea can
   still be a good one — what has to change is either `opts.bound.flag`,
   which turns the bounding off, or the definition, so that the middle half of
   the column is not one value.
2. **Correlation with the performance values.** This is the check SIFTED
   itself applies, so it decides whether the feature survives selection at
   all — which means reproducing what SIFTED does rather than improving on
   it. Pearson, on the *processed* features, against `model.data.Y`: PRELIM
   runs first and overwrites `data.X` with the bounded, Box-Cox'd, z-scored
   values, and SIFTED correlates those. `Y` is raw under an absolute
   threshold, the deficit from the best algorithm under a relative one. An
   algorithm that is good on no instance is not in `Y` at all, because the
   build removes it before PRELIM normalises (`InstanceSpace.m:636-645`), so
   a correlation with it counts for nothing; the script leaves such
   algorithms out and names them. Keep
   anything reaching |rho| >= `opts.sifted.rho` at p <= `opts.sifted.pval` for
   any one algorithm. Those default to 0.1 and 0.05, but where an existing
   space is in play read its own `model.opts.sifted` instead of assuming — the
   run may have set them differently. Clearing the floor is sufficient but not
   necessary: SIFTED keeps whichever feature correlates most with each
   algorithm whatever that correlation is (`core/SIFTED.m:121`), so a candidate
   under the floor can still survive by being the best there is.

   Before discarding anything here, check whether the relationship is
   non-monotone. Pearson scores a U-shape near zero, so a feature that
   predicts performance well on both sides of a turning point fails this check
   at a correlation near nothing. The script fits a quadratic to everything the
   check drops and reports any that reach the same floor. The remedy is to
   redefine the feature so the relationship runs one way — the distance from
   the turning point rather than the raw value — because the pipeline is
   Pearson from end to end and will not find it otherwise.

## Rank

Everything that survives gets scored and ordered. The measurements below
answer different questions and can disagree.

**Both correlations are computed per algorithm**, one number for each. A
candidate has to clear the floor on one algorithm rather than on all of them,
and which algorithms it explains is part of what you are reporting, so keep
the whole vector rather than only a summary of it. Ordering the candidates
does need a single number, and the largest |r| across the algorithms does that
job — name the algorithm it came from when you use it.

- **Correlation with the performance values**, already computed by the check
  above, one per algorithm. This is what decides whether SIFTED keeps the
  feature, so it comes first among equals, and larger magnitudes are what you
  want.
- **Correlation with the binary good/bad label** — point-biserial against
  `Ybin`, again one per algorithm. This is what PYTHIA predicts and what TRACE
  draws footprints from, so a candidate can clear the correlation check and
  still be the weaker of two for everything the pipeline does after selection.
- **Spread across the portfolio** — the distance from a feature's strongest
  correlation to its weakest. Wide spread means the feature separates the
  algorithms; near-zero spread means it tracks difficulty, which tells you
  what is hard rather than which algorithm to pick. Prefer wide, and check
  that the ranked set is not uniformly narrow — though this one counts for a
  bit less than the others.

Two more where a feature set already exists. They are separate measurements
answering separate questions, so report both rather than running them together:

- **Overlap with what is already there** — max |r| against every existing
  column, naming the nearest match. This says whether you have written
  something the set already measures. A duplicate does no harm, since SIFTED
  clusters the two together and keeps one of them, but it buys nothing either,
  so it ranks last whatever its own correlation. The script reports the
  clusters as well, so which column a candidate is competing with can be read
  off rather than inferred from a correlation.
- **Incremental value** — the cross-validated predictive gain from adding the
  candidate to the existing set. Overlap says whether the measurement is new;
  this says whether the newness is worth anything, and it says it out of
  sample. Read the gain against the headroom rather than in absolute terms,
  because on a set already scoring near the ceiling a small gain is a ceiling
  effect rather than evidence the candidate is worthless.

`scripts/evaluate.py` runs both checks and every measurement in this section.

## Output

One record per proposal, and no bare numbers:

- **Name**, and a plain-language definition — one line where the feature is
  simple enough for that, longer where it is not.
- **Code** against the documented instance format.
- **Why it should matter** — what the feature captures, and why that affects
  this algorithm. With a citation where one exists.
- **Measured sign**, per algorithm where they differ. Mark `D` for a feature
  that *differentiates algorithms* rather than tracking difficulty.
- **What it was written for** — the audit finding, the region nothing
  explains, or the algorithm nothing currently accounts for. In a blind pass
  there may be no such target, and saying so beats inventing one.
- **Self-test**: runs on every instance without erroring, and is not constant
  across them.

Rejects get a record too, with the reason. ISA is iterative, and a feature that
fails now can pass once the instance set covers more of the space.

Then the ranked table: one row per survivor, carrying its correlation with
the performance values, its correlation with the binary label, the algorithm
each of those belongs to, and its spread across the portfolio. Where a feature
set already exists, the overlap and the incremental gain go beside it as their
own tables — the gain is one figure per algorithm rather than per feature, so
it cannot be a column here. Say which of them you ordered on, because they can
disagree and the reader needs to know which one you let decide.

