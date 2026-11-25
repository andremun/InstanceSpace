---
name: instance-space-improvement
description: >
  Decide what to change about an instance space analysis: which empty
  regions are worth filling and what the targets in them should be, and
  whether the feature set, the algorithm portfolio, the performance
  measure or the pipeline options is the real limitation. Use when the
  user has a built space and wants recommendations rather than a reading
  of it. Triggers on "how do I improve my instance space", "where are
  the gaps in my benchmark", "what instances should I generate", "should
  I add another algorithm", "is my feature set the problem", "what is my
  analysis missing", "what should I do next with this space". For
  reading the space in the first place use
  `instance-space-interpretation`, which this skill assumes has been
  done. For producing instances at the targets this skill chooses, use
  `instance-generation`.
---

# Improving an instance space

The space is built and read. The job now is to say what to change, and
which change is worth making first.

<!-- Boundary with instance-generation: the screen here is eliminate-only,
     and nothing downstream tests reachability -- a search scored on
     closeness to the target does that work instead. What a target looks
     like in raw units belongs to that skill and is not to be restated
     here. -->

Four adjacent skills own the remedies rather than the diagnosis:
`feature-proposal` works out what a missing feature should be,
`performance-thresholds` settles the good/bad rule, `instance-generation`
produces instances at targets this skill hands it, and
`instance-space-analysis` covers the pipeline options. This skill decides
which of those is the live problem.

A space can be inadequate in five ways — its instances, its features,
its portfolio, its performance measure, its pipeline options. The
evidence for each comes from Steps 1-8 of
`instance-space-interpretation`; this skill collects it and says which is
worth acting on.

## Instances: the empty regions

**Check what was held constant across every instance.** Read the
instance files or the collection's documentation, not only the feature
table: a variable fixed at one value has no gradient, so the projection
cannot show it and no amount of staring at the map will reveal it. The
candidates are problem size, time or resource budgets, and any generator
parameter varied in one dimension but not another. If a size-like
variable is fixed, say so plainly — the study cannot speak to how
performance scales, which is usually a reader's first question. Then
check whether instances already exist that were never used: an archive
with unused instances is the cheapest fill there is, needing solver runs
rather than generation.

**Measure the coverage.** Compare the instances against the outer
boundary, preferring the correlation-pruned one — `model.cloist.Zecorr`
or `bounds_prunned.csv`, which keeps only the corners consistent with
the observed pairwise feature correlations. `model.cloist.Zedge` and
`bounds.csv` use the feature ranges alone and are wider. Report the hull
around the instances as a share of the boundary's area, and the count of
empty cells on a grid restricted to cells inside the boundary — the
second finds interior holes the first hides. Check first which boundary
you have, as Step 1 of `instance-space-interpretation` describes: above
`opts.cloister.maxFeatures`, 20 by default, both are the hull of the
instances themselves and every coverage number becomes circular
(`core/CLOISTER.m:80-88`).

**Rank what you would fill.** Highest value is a gap beside a region
where one algorithm wins alone, since filling it tests the report's
central claim; then a gap beside a region carried by a single instance
source, which needs independent coverage before it can be trusted;
lowest, a gap inside a region where everything works or nothing does,
since neither separates the algorithms. For each, report how the
algorithms behave on the nearest existing instances — that is a
falsifiable prediction, and the first thing to check if the instances get
made. A list of empty coordinates is not a recommendation.

**Screen the gaps against the boundary, in the one direction it is good
for.** A gap outside the boundary needs some feature to take a value
beyond anything yet observed, which is enough to drop it from the
ranking. A gap inside tells you nothing: the boundary comes from the
feature ranges and the *pairwise* correlations, so it routinely marks
impossible regions reachable. The screen eliminates and never admits,
and above the feature cap it is not worth running: the boundary is then
the hull of the instances, so it rejects every gap outside the coverage
you already have and tells you only that the gap is empty. Every
surviving gap is therefore still unresolved, and nothing downstream will
resolve it either: `instance-generation` aims a search at the target and
reports how close it got, which is the answer you get. So rank a gap on
what filling it would settle, say that the screen did not rule it out,
and let the miss distance speak afterwards.

**Then say what the targets are, not just which regions.**
`instance-generation` takes targets as given and will not choose them,
so a recommendation that stops at a region is unfinished. Four things to
settle.

- **Which space each target lives in.** A point in the two-dimensional
  space leaves the remaining feature directions free: good for covering
  the map, no control over the features themselves. A full feature
  vector pins every feature — tighter control, harder search. Neither is
  wrong, but which one you mean has to be said: the two are not
  comparable once the misses are reported.
- **Where to sample them.** Sampling inside the feature ranges you
  already have mostly returns the coverage you already have. Munoz &
  Smith-Miles (2020) found a Latin hypercube over the observed feature
  space produced instances that "do not push the boundaries of the
  existing benchmarks very far", while sampling the two-dimensional
  space over bounds deliberately wider than the observed ones is what
  reached the edges. If the aim is to extend the space, aim past the
  current hull — and expect those instances to need a rebuilt space
  before they can be seen where they really are, since the stored model
  clips anything beyond the feature ranges it was fitted on.
- **Several, not one.** A single point buys one instance and no way to
  tell a reachable neighbourhood from a lucky hit. Spread them across
  the region, far enough apart to be distinguishable in the projection.
- **Location or behaviour.** Targets at locations fill a space.
  Targets defined by algorithm behaviour — instances evolved to be hard
  or easy for one algorithm (Smith-Miles et al., 2010), or to maximise
  the gap between two — probe what an algorithm is sensitive to instead,
  and Munoz & Smith-Miles (2020) found such instances concentrated in
  one part of the space rather than spread across it. They also cost far
  more, because scoring a candidate needs an algorithm run rather than a
  feature computation. Use them to stress-test a claim about an
  algorithm, and location targets for coverage.

## Features

Step 4 listed what selection dropped, Step 5 sorted the survivors into
clear and muddy, and Step 7 may have found a region the features do not
explain at all. All three are feature findings rather than instance
findings, and each one needs a different fix, so keep them apart before
you recommend anything.

The evidence sits in `model.sifted.rho` and `model.sifted.p`, the
feature-by-algorithm correlation and p-value matrices, one row per
feature. Both are computed over every original feature, before anything
was dropped (`core/SIFTED.m:115`), so you can read them for the features
that did not survive as well as the ones that did. The two names are
close and mean different things: `opts.sifted.rho` is the threshold,
`model.sifted.rho` the matrix.

- **A survivor the projection could not show** is the strongest finding
  here. The projection is a compromise across every surviving feature,
  so a muddy plot means the picture cannot show that feature, not that
  the feature does not matter. One that correlates with performance and
  still cannot be drawn is a reason to revisit the set, and not a reason
  to generate instances.
- **A feature the correlation filter dropped** was weakly correlated
  with every algorithm — that is what the filter tests — so there is
  usually nothing to fix. `feature-proposal` has the rule itself and the
  one case worth re-checking. `model.sifted.p` shows a second: a real
  effect that missed the significance cut because the instance set was
  small.
- **A feature the clustering stage dropped** cleared the filter first,
  then lost to another feature in its own cluster. The search keeps one
  feature per cluster and scores them on how well they predict the
  good/bad labels, by k-nearest-neighbour accuracy
  (`core/SIFTED.m:147-203`), so the survivor predicts at
  least as well. That is redundancy, not irrelevance: the set is bigger
  than it needs to be rather than missing something. To tell this case
  from the one above, read the dropped feature's row of
  `model.sifted.rho` against the run's own `opts.sifted.rho`.

`model.sifted.clust` holds the cluster membership, so you can name the
survivor standing in for a dropped feature — but its rows are the
features that passed the filter, and `selvars` was overwritten with the
final survivors afterwards (`core/SIFTED.m:204`). Re-apply the filter
rule to `model.sifted.rho` to recover that row order.

Where a region's performance has no feature story at all, name the
quantity you would need to measure instead. `feature-proposal` is the
skill for working out what that feature should be.

## Algorithms

Step 6's regions and Step 8's per-algorithm models say whether the
portfolio is doing any work. `model.trace.summary` has the numbers,
giving each algorithm two footprints as a fraction of the whole space.

- **Dropping one.** A small `Area_Good_Normalized` says the algorithm is
  rarely good; a small `Area_Best_Normalized` beside a healthy good
  footprint says it performs well but is never the best choice, which is
  redundancy rather than weakness. Both columns take their prefix from
  the measure, so in a three-dimensional space they are `Volume_`
  instead (`core/TRACE.m:107`). Before you treat either as zero, check
  how many instances the footprint rests on. The rule depends on the
  TRACE method. Under TRACE3, a footprint resting on no more than
  `opts.trace.minInstances` instances, 4 by default and counted over the
  unique instances the algorithm is both predicted good and actually
  good on, or covering less than `minAreaFrac` of the space, is thrown
  out with a measure of 0 (`core/TRACE.m:265-281`). The legacy method
  has neither setting. It drops a footprint only when its instances sit
  at fewer than three distinct points (`core/TRACE_legacy.m:115`).
  Either way, a dropped footprint looks exactly like having no region.
  If the evidence is instead that the classifier predicts good almost
  nowhere, check `model.pythia.recall` first: a classifier that
  collapsed onto the majority class does that whatever the algorithm
  does.
- **Adding one.** A region where every algorithm fails is the case for a
  new algorithm, and those instances are `model.data.numGoodAlgos == 0`.
  `model.trace.hard` is not a substitute: it is the footprint of
  instances where at most `opts.perf.betaThreshold` of the portfolio is
  good, 55% by default (`core/PRELIM.m:145`). Only an absolute threshold
  can show this at all — a relative one makes the best algorithm good on
  every instance by construction — so check which was used before
  concluding there is no gap.
- **Redundancy.** Two unrelated methods good in the same region means
  the portfolio is carrying duplication there. Read it off the *good*
  footprints, not the best ones. Whether best
  footprints overlap depends on which TRACE ran: the legacy one removes
  the overlap between them, and TRACE3, the default, leaves it
  (`core/TRACE.m:176-178`), so under TRACE3 that overlap can be an
  artefact of the method rather than a finding. Read
  `model.opts.trace.method` rather than assuming the default.

## The performance measure

The good/bad threshold decides every region in the picture, so start by
naming the one that was used: `model.opts.perf` holds `AbsPerf`,
`MaxPerf` and `epsilon`. On a model saved before `MaxPerf` existed you
may find only the legacy `MaxMin`, or neither
(`InstanceSpace.m:968-975`).

Then work out whether it is the live problem. A good-rate across
`model.data.Ybin` (`algorithm_bin.csv`) near 0 or 1 means the labels
carry almost nothing; a per-algorithm rate of exactly 1 means that
algorithm contributed nothing to the picture. An algorithm good on no
instance does not show up as a rate of 0, because the build removed it
before anything else ran (`InstanceSpace.m:636-645`). Compare
`model.data.algolabels` with the metadata's `algo_` columns to find
one. Either way the threshold
is the first thing to change, not the instances. A censored runtime is
the other common case, mixing slow with never finished, so the label
says less than it appears to.

Say what the interpretation would look like under a different threshold.
`performance-thresholds` is the skill for choosing and justifying one.

## The pipeline options

Some settings limit what the analysis can say at all.

- **Feature selection may not have run.** With three or fewer features
  on input, SIFTED skips both of its stages and keeps everything
  (`core/SIFTED.m:105-112`). After the correlation filter it skips the
  clustering and search stage whenever three or fewer features survive,
  and also whenever the survivors do not outnumber `opts.sifted.K`, 10
  by default (`core/SIFTED.m:132-145`). The second is the one that
  usually bites, and it means the feature set is whatever the filter
  left.
- **The feature clustering may have been poor.**
  `model.sifted.Ksuggested` names a better cluster count, and is written
  only when the silhouette for the K that ran was bad enough to warrant
  one, with `opts.sifted.diagnostics` on as it is by default
  (`core/SIFTED.m:156`). Check the gate above before reading anything
  into its absence: where that skipped the stage there was no clustering
  to judge.

Where a setting has constrained the conclusions, say which and what it
would take to lift it. `instance-space-analysis` covers the options
themselves.

## Then say which one you would do first

And what it would settle. Every dimension above is cheaper than
generating instances, so a report that recommends instance generation
without having weighed them is recommending the most expensive option by
default.
