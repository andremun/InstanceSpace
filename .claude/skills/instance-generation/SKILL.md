---
name: instance-generation
description: >
  Generate new problem instances that land at target locations already
  chosen in an existing instance space: steer a generator or a search
  towards a target, and check where the instances actually landed. Use
  when the targets are settled — by `instance-space-improvement` or by the
  user — and the question is how to reach them. Also covers labelling the
  new instances and folding them back into the analysis. Triggers on
  "generate instances", "how do I make instances for this gap", "fill the
  empty regions", "target instances at this point", "evolve new
  instances", "what instance would land here", "new benchmark instances",
  "instance generator". For choosing the targets, deciding which gaps are
  worth filling, and whether instances are the right thing to change at
  all, use `instance-space-improvement` first.
---

# Generating instances at a target location

The space is built, it has been read, and targets have been chosen. This
skill turns them into instances.

`instance-space-improvement` decides which gaps matter, whether
instances are the right lever at all, and what the targets should be —
run it first, because everything else you could change is cheaper than
building a benchmark. **This skill does not choose targets.** It takes
them as given, from that skill or from the user, and owns the rest: how
to reach a target, and how to tell whether it arrived. New feature
definitions come from `feature-proposal` and existing ones whose code
has gone missing from `feature-recovery`, the good/bad rule from
`performance-thresholds`, and the pipeline from
`instance-space-analysis`.

## Inputs

Ask for anything you are missing before you start. Send one message,
wait for the reply, then get on with it using whatever you have.

**Required**

- **The built model** — `model.mat`, or the output CSVs together with the
  projection matrix. The projection, the stored preprocessing
  parameters and the boundary all come from here.
- **The targets, and what they specify.** From
  `instance-space-improvement`'s ranked gaps, or from the user. Two
  things have to come with them, since this skill does not decide them:
  what each target asks for, and how many there are. A target is either
  a location — a point in the two-dimensional space, or a full feature
  vector — or a behaviour, such as an instance one algorithm finds hard,
  or one that separates two algorithms. Ask, if a target arrives as a
  bare pair of coordinates. With `k` features and two axes, a 2-D target
  leaves `k-2` feature directions free, so the kinds are scored
  differently in §2 and their misses cannot be compared in §3.
- **The instance format, and what makes an instance legal** — the
  constraints a solver will reject a candidate for. §1 chooses whether
  to search over a generator's parameters or over the instances, and the
  generator or the operator can be written as part of the work; what
  cannot be worked out from the model is what a valid instance looks
  like. Example instances carry the format implicitly, so it is usually
  the legality rules that have to be asked for.
- **Code that computes every feature** on a new instance, matching the
  definitions behind the model. The projection takes a full feature
  vector, so a candidate you cannot measure cannot be scored. Where the
  code is gone but the instances and their feature table are not,
  `feature-recovery` reconstructs the definitions from the values — and
  its own standard, that a definition reproduces the existing column, is
  the one to insist on before running it on anything new.
- **The existing instance files**, not only their coordinates. A
  perturbation operator needs seeds, §2 may start from the instances
  nearest the target, §1 wants any that were never used, and the
  recovery route above reads them. The model carries the feature table
  but not the instances themselves.

**Optional**

- **An existing generator**, with the parameters you can set. §1 can
  search over the representation without one, so this is a saving rather
  than a requirement — and the saving is not only effort: a published
  generator is code nobody has to take on trust, which matters for
  instances other people are meant to use.
- **The algorithms and the run harness**, with the budget and the
  hardware the original runs used. Needed to label the new instances.
  Without it you can still deliver instances at the target, but not
  their performance.
- **The reference value the performance measure needs** — an optimum,
  best known solution or bound, per instance. Usually missing for
  anything newly generated; §4 says what to do about that.
- **Which features are set at generation time and which are only
  measured afterwards.** The one most worth asking for: it decides
  whether §1 can reach a target by setting a parameter rather than
  searching for one.

### When something is missing

Ask for all of it at once, say in a few words what each thing would let
you do, and make clear what you will carry on without. Then stop and
wait.

Two of the required inputs cannot be worked around. Without the model
there is no projection, no stored preprocessing parameters and no
instance coordinates, so nothing here runs at all. Without targets there
is nothing to aim at, and this skill does not choose them —
`instance-space-improvement` does. If either is still missing after you
have asked, say that is why you stopped, rather than generating
something unaimed.

For the rest, one round is enough. Having no generator does not stop
anything, since §1 can search over the representation instead, but not
knowing what makes an instance legal does: an operator free to produce
candidates the solvers reject spends the whole budget on nothing.
Missing feature code can be rebuilt wherever the instances and their
feature table survive, at the cost of a recovery pass first; say that it
cost one. Missing the instance files closes that route, leaves a
perturbation operator with no seeds and §2 with nothing to start from.
In any of those cases say what is missing and what it costs rather than
delivering instances nobody can place: without feature code in
particular there is no way to report where anything landed, which is
most of the result.

## 1. Choose how to produce the instances

Two things you can search over: the generator that makes instances, or
the instances themselves. The choice is between those two, never between
generating and not generating — `instance-space-improvement` has already
weighed an existing archive against generating and come down on
generating, so that question is settled before this skill starts.

**Adjust or evolve a generator.** Where one exists, or where you can
write one, search over its parameters rather than over instances. The
search space is far smaller, each evaluation is a generator run plus a
feature computation, and what comes out is a recipe rather than a single
file — re-run it on other seeds for many instances at the same target,
which is what a benchmark wants. Sometimes there is no search to do at
all: where one feature climbs steadily as one parameter rises, just try
values of that parameter until the feature reaches the target. Its limit
is reach: you get only the instances the parameterisation can express,
so a target needing structure the generator cannot produce stays out of
reach however well you tune it. Parameters that control structure are
what make it worth doing: benchmark libraries and randomly generated
instances tend to sit together in the middle of the space at average
feature values, while instances from generators with structural
parameters spread across it.

**Generate or evolve individual instances.** This is for the features no
parameterisation reaches, which is the usual case: there is often no way
to know how to adjust a random instance generator so that it produces
instances with a particular value of a measurable feature. Search over
the instances themselves, with whatever optimiser suits the
representation. A genetic algorithm over the instance format is the
usual choice, or genetic programming where the instance *is* a function,
but nothing in the method depends on either: simulated annealing, tabu
or plain local search need only a neighbour operator, which is the
easier fit for a representation that is awkward to cross over, and a
continuous parameterisation suits CMA-ES or Bayesian optimisation,
worth it when one feature evaluation is expensive. Choose on the
representation and the cost of an evaluation, then say which you used.
§2 covers running either search.

Whatever the operator, it must return a legal instance: a perturbation
that produces something the solvers reject is a wasted evaluation, and a
run of them looks like a search that cannot reach the target. Work out
what legal means from the problem's definition rather than from the
feature table — a relation every existing instance satisfies may be a
limit of the one generator that made them rather than of the problem,
and testing a candidate constraint against the real instances only
catches the ones they contradict. The relations to trust are arithmetic
ones between features measuring the same quantity, where a minimum
cannot exceed a median, and shares of one whole, which cannot sum past
it.

**Worth doing alongside either.** Compute the features of any instances
that were never used, project them, and see where they fall. It is cheap
and it pays twice — instances near the target are seeds §2 can start
from, and instances sitting *in* the target region mean the region was
never as empty as the map suggested, which belongs in the report.
Neither finding replaces generating.

**Turn the target into feature values, if the method needs them.** Any
search scored on distance in the two-dimensional space does without
this, whether it searches instances or parameters. Setting a parameter
to hit a feature value needs it, because there you are aiming one
feature at a number rather than measuring a distance. The projection is
linear, so the smallest of the feature vectors that project onto a
target `z` is `x = Xmean + A' (A A')^-1 z`, in processed units. `Xmean`
is `model.pilot.Xmean`, present only for the analytic projection and
zeros otherwise — the forward projection is `Z = (X - Xmean)*A'`
(`InstanceSpace.m:1036`), and dropping the shift puts the target in the
wrong place on a model that has one. Verify the projection reproduces
the stored coordinates before inverting anything; if it does not, the
matrix and the feature table are in different feature orders. Then
invert PRELIM one feature at a time, because nobody can build an
instance from a z-score: un-z-score with `muX`/`sigmaX`, invert Box-Cox
at `lambdaX` (`y = (x^L - 1)/L`, so `x = (L*y + 1)^(1/L)`, and `exp(y)`
at `L = 0`), subtract 1 and add back `minX` (`core/PRELIM.m:241-248`).
The parameters are indexed over every feature PRELIM saw, so pick out
the selected ones with `model.featsel.idx`.
With `k` features and two axes a point leaves `k-2` directions free, so
say which features the value you picked actually pins, and which of them
can be controlled at generation time.

## 2. Run the search

**The objective is distance to the target.** Minimise `||v - v_T||`.
You will also see it written as maximising `1/(||v - v_T|| + e)` for a
small `e`, which is the same problem turned round to give a genetic
algorithm something to call fitness; with any other optimiser, use the
distance.

**A behaviour target scores differently, and costs far more.** Where the
target is an instance one algorithm finds hard, or one that separates
two algorithms, the objective is the measured behaviour itself — the
algorithm's result on a candidate, or the gap between two results —
rather than a distance in the space. The rest of this section still
applies, but each evaluation now needs an algorithm run instead of a
feature computation, so expect far fewer of them: set the budget before
you start, and keep the starting points few and well chosen. Measure
each kept candidate's location as well, by the rule below, since that is
what §3 reports.

**Compute `v` the way the model did.** The projection acts on
*processed* features. A candidate's raw features must go through the
stored transform — clip to `[lobound, hibound]`, shift by `minX`,
Box-Cox at `lambdaX`, z-score with `muX`/`sigmaX` — never a transform
refitted on the candidates, which would put `v` in different units from
the model. The toolkit does this itself: PRELIM takes a fourth argument
and switches to evaluation mode, applying a stored run's parameters
instead of fitting new ones (`core/PRELIM.m:63,151`). Subset the
features with `model.featsel.idx` afterwards, which means the full
feature set has to be computed even though only the selected ones are
projected.

**Pick the starting points on purpose, and use more than one** — a
population or a single state, depending on the optimiser. A random start
keeps diversity but begins inside the cloud of instances you already
have. Starting from the existing instances nearest the target converges
faster at some cost in diversity — worth it when the target is far out,
and worth avoiding when the point of the exercise is instances unlike
the ones you have. Mixing them works well: half the population drawn at
random, half the existing instances closest to the target.

A population carries several starts already; a single-state method does
not, so restart it from different points and keep the best. A search
over instances has many local optima, so one run that never approached
the target tells you more about the run than about the target — and §3
will ask you to tell those apart.

**Stop on the miss distance, not the iteration count.** Set a threshold
on it and stop there. And keep every candidate that came close, not just
whatever the search finished on. A converged population is a cluster of
near-identical instances; a single-state method throws away everything
it passed through on the way. Take only the last answer and the whole
budget has bought one point of coverage.

## 3. Check where the instances landed

**Report the miss distance, aggregated and per target.** The usual
summary is the median RMSE per feature across a batch of targets, with
its interquartile range. Say which space you measured it in, since a 2-D
miss and a feature-vector miss are different quantities. Plot the
targets and the achieved locations together. If the misses all pull the
same way that is worth reporting.

**For a behaviour target, report the behaviour and the location.** The
behaviour you asked for against the one you got — how hard the instances
actually were, or how large the gap turned out — and where they landed,
since asking for a behaviour says nothing about where in the space the
instances will sit. Instances evolved for behaviour tend to cluster in
one part of the space rather than spread across it, so report the spread
rather than assuming it.

**A target the search never reached is not a target nothing can reach.**
Failing to produce an instance at a point means the search could not
find one, not that none exists — it may need tuning or more attempts.
Report a miss as a miss, with the budget it was given, and leave the
impossibility claim alone.

**Check the new instances are distinct.** From each other, and from the
instances already in the space. Near-duplicates inflate the count
without adding coverage, and they bias every later feature correlation
towards whatever the duplicated instance happens to be.

## 4. Fold them back in

**Run the algorithms under the original conditions** — same budget, same
hardware, same measure — but only where the runs are definitely cheap,
such as a constructive heuristic. Otherwise do not measure performance
at all: an unlabelled instance still carries its location, which is the
result this skill set out to produce. Never economise on the conditions
themselves, because an instance measured differently is worse than one
left unmeasured — the difference will read as a real one. A relative
measure is the common reason you cannot label at all: it needs an
optimum, a best known solution or a bound, and a freshly generated
instance has none of them.

**Project, then rebuild, and keep the two apart.** Pushing the new
instances through the stored model keeps them comparable with the
original coordinates, which is what §3's check needs; the toolkit's
`explore` does this from a `metadata_test.csv` in a test directory,
reusing the stored preprocessing, projection, selector and footprints
and leaving the boundary alone (`InstanceSpace.m:952`). Its `feature_`
columns have to be the same ones as the training metadata and in the
same order: INIT checks the count against the fitted transform width and
the order against `model.featsel.labels`, and errors on either
(`core/INIT.m:138-166`). On a model old enough to lack that label list
the order check is skipped, and columns are silently clipped and
normalised with the wrong feature's parameters instead. Rebuilding
instead — rerunning the pipeline with the new instances included —
refits the preprocessing, may select different features and produces a
different projection, so the old coordinates no longer apply. Do the
first to verify the instances arrived, the second to get the improved
analysis, and never put both in one figure.

**Watch for the clipping warning.** Evaluation mode warns when more than
5% of the new instances get clipped to the training feature bounds
(`core/PRELIM.m:166-171`). Expect it when you have aimed at empty
regions, and read it rather than dismissing it: clipping pulls those
instances back inside the old distribution, so no search however good
can reach a target that lies beyond those bounds once projected.

**Check the prediction, where you labelled them.**
`instance-space-improvement` ranks a gap with a prediction of how the
algorithms behave on the nearest existing instances, and a labelled
instance at the target is what tests it. One that performs as predicted
confirms the region. One that performs differently tells you more: the
features are missing whatever actually matters there, and fixing that is
a `feature-proposal` job rather than a generation one.

## Reporting

The report needs all five:

- the targets, and why those ones;
- how the instances were produced — which method from §1, and the
  optimiser or the parameter you set;
- how many were made, and how close each came to its target;
- whether they were labelled and under what conditions, or that they
  were not;
- what filling the region settled, if they are labelled.

Instance files on their own, with no miss distances, are not something
anyone else can build on.
