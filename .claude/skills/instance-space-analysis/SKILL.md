---
name: instance-space-analysis
description: >
  Perform, interpret, or advise on Instance Space Analysis (ISA) -- the
  methodology for objective algorithm benchmarking via a 2D/3D projection
  of problem instances (Smith-Miles & Munoz, ACM Comput. Surv. 55(12),
  2023). Use this skill whenever the user wants to: build an instance
  space from feature/performance meta-data; interpret PRELIM, SIFTED,
  PILOT, CLOISTER, PYTHIA, or TRACE/TRACE3 output; design a new ISA
  application for a problem domain; work with this repository's MATLAB
  `InstanceSpace` toolkit (build/explore, options.json, the InstanceSpace
  class); or assess whether a benchmark suite or instance space is
  adequate. Trigger on "instance space", "algorithm footprint", "PILOT
  projection", "TRACE3", "MATILDA", "ISA toolkit", or a request to
  run/interpret this pipeline on new meta-data. See
  references/matlab-toolkit.md for operational (file-format, class API,
  options.json schema) detail specific to this repository; this file
  covers methodology, decisions, and pitfalls.
---

# Instance Space Analysis (ISA)

ISA extends Rice's Algorithm Selection Problem (1976) with a geometric
layer: instead of learning a selection mapping directly from features to
algorithms, it projects instances into a 2D (or 3D) plane and identifies
contiguous regions -- **footprints** -- where each algorithm performs
well. The result is not just a predictor; it is a visual, falsifiable
account of *why* an algorithm wins where it wins, and whether the
benchmark suite used to draw that conclusion was diverse enough to trust.

Primary source: Smith-Miles & Munoz, "Instance Space Analysis for
Algorithm Testing: Methodology and Software Tools," *ACM Computing
Surveys* 55(12), Article 255, 2023 (the tutorial this skill is built
from). Supersedes/extends: Munoz, Villanova, Baatar & Smith-Miles,
*Machine Learning* 107(1), 2018 (PILOT); Simpson, Munoz, Kandanaarachchi
& Campello, "ISA3: a 3-dimensional expansion of instance space analysis,"
*Machine Learning* 114:240, 2025 (3D projection, TRACE3).

This skill is scoped to the `andremun/InstanceSpace` MATLAB toolkit --
the reference implementation and the engine behind MATILDA. Everything in
this file that describes *current code behaviour* was verified directly
against that repository (v0.9.0) at the time this skill was written, not
reconstructed from the papers. Where the shipped code has moved past what
a paper describes, that is flagged explicitly -- treat the code as ground
truth per the source-of-truth discipline in the persona/
scientific-computing skills, and re-verify against the live repository
before trusting a detail below if it might have changed since (this
toolkit is under active refactor; option names and defaults have already
changed multiple times across recent versions -- see the repo's
`RELEASE_NOTES.md` for the version you are actually looking at).

---

## 1. The conceptual framework -- six spaces

| Space | Symbol | Meaning |
|---|---|---|
| Problem space | P | All relevant instances of the problem, most of which are never observed |
| Instance subset | I ⊂ P | The instances actually collected, with meta-data |
| Feature space | F | Each instance x is a vector f(x) characterising its structural difficulty |
| Algorithm space | A | The portfolio of algorithms compared |
| Performance space | Y | A user-defined "goodness" measure y(alpha, x) per (algorithm, instance) pair |
| Instance space | Z ⊂ R^2 (or R^3) | The projected 2D/3D coordinates the whole analysis is read from |

Meta-data is two matrices: F ∈ R^(m×n) (m features, n instances) and
Y ∈ R^(a×n) (a algorithms). Rice's original framework only learns
S(f(x), y) → alpha*; ISA additionally learns a projection g(f(x), y) → z,
so that *why* is answered geometrically, not just *which*.

---

## 2. The six-step iterative methodology

1. **Collect meta-data**: instances I, features F, performance Y for a
   portfolio A. This step is the one most likely to determine whether
   later insights are trustworthy -- see Section 6.
2. **Construct the instance space**: feature selection, projection, and
   theoretical boundary (PRELIM -> SIFTED -> PILOT -> CLOISTER; Section 3).
3. **Automated algorithm selection**: learn alpha* from instance
   coordinates (PYTHIA; Section 3).
4. **Generate footprints and metrics** (TRACE / TRACE3; Section 3).
5. **Analyse**: distribution of sources, classifier accuracy, footprints,
   feature distributions (Section 4).
6. **Augment**: generate new instances, algorithms, or features to close
   gaps identified in Step 5; return to Step 2, or stop if the space has
   converged (no new features selected, no gaps in the boundary).

Convergence is a real stopping criterion, not a formality: ISA is
finished for a given algorithm set when the boundary's interior is
densely filled by real or generated instances and the feature set is
stable across an added iteration. Treat "we ran ISA once" as an
interim result, not a final one, unless the domain genuinely has no route
to generate more instances.

---

## 3. The five core methods

### PRELIM -- Preparation for Learning of Instance Meta-Data

Two jobs: (a) define binary "good performance," and (b) bound and
normalise the meta-data.

- **Good performance** is user-defined and consequential, not a detail:
  *absolute* (`opts.perf.AbsPerf = true`: performance passes a fixed
  threshold `opts.perf.epsilon`) or *relative* (within a margin of the
  best algorithm on that instance). `opts.perf.MaxPerf` separately
  selects whether performance is maximised (efficiency) or minimised
  (cost). The paper calls the absolute-vs-relative choice "somewhat
  arbitrary but exerting significant influence on the final results" --
  do not let it be an unexamined default; state and justify it the way
  you would a p-value threshold.
- Bounding: each feature clipped to median +/- `opts.prelim.iqrMultiplier`
  (default 5) times the IQR (robust to outliers by construction). A
  feature with more than `opts.prelim.nanThreshold` (default 0.20)
  fraction of missing values is dropped entirely before bounding.
- Normalisation: one-parameter Box-Cox (lambda fit by maximum
  log-likelihood) then z-score, applied to both features and
  performance.
- **Code detail beyond the paper**: `PRELIM.m` also converts performance
  into a *relative deficit* from the best algorithm on each instance
  (`Y = 1 - Y/Ybest` maximising, `Y/Ybest - 1` minimising) before
  normalisation -- not the raw performance values. It also derives a
  per-instance `beta` (easy/hard) flag from `opts.perf.betaThreshold`
  (default 0.55) fraction of algorithms achieving good performance, used
  later by TRACE's beta-hard footprint. Ties in "best algorithm" are
  broken at random per instance, logged as a percentage.

### SIFTED -- Selection of Instance Features to Explain Difficulty

Goal: a small, uncorrelated feature subset that explains algorithm
performance. Implemented as `SIFTED.m` (the historical `SIFTED2.m` name
now lives under `deprecated/` as a warn-and-forward alias).

- Textbook version (2023 tutorial, Algorithm 2): correlation filter
  (keep the top feature per algorithm, plus any feature with |correlation|
  >= 0.3) -> k-means clustering on 1-|correlation| dissimilarity (K
  clusters, default 10) -> one feature per cluster chosen by lowest
  out-of-bag error of a Random Forest on a temporary PCA 2D projection.
- **Code has moved past this description.** `SIFTED.m` still does the
  correlation filter (`opts.sifted.rho` default 0.1 minimum correlation,
  `opts.sifted.pval` default 0.05 significance -- not the paper's 0.3
  threshold), and still clusters on 1-|correlation| with k-means
  (`opts.sifted.K` default 10, `opts.sifted.MaxIter`/`.Replicates`
  controlling convergence), but replaces the PCA+RandomForest
  representative-selection step with a **genetic algorithm** whose
  fitness is the worst-case (maximum, across algorithms) k-fold
  cross-validated KNN classification loss for a candidate one-feature-
  per-cluster combination, projected via **PILOT's analytic branch** (a
  single closed-form solve, no BFGS restarts) at the same dimensionality
  as the outer pipeline's `opts.pilot.dims`. Previously-evaluated
  combinations are cached (a `persistent` variable) for the duration of
  the call. This is a materially different (and more expensive) search
  than the paper describes. If reproducing a published number or a
  paper's exact SIFTED description, verify which version and which
  commit produced it.
- Edge cases handled explicitly in code: <=1 feature is a hard error;
  <=3 features skips the clustering+GA stage entirely (correlation-filter
  output is the final answer); fewer features than K skips clustering. A
  silhouette-score advisory suggests a better K if the chosen K clusters
  poorly (below 0.5) -- read this warning, since a bad K silently
  degrades the whole downstream projection.

### PILOT -- Projecting Instances with Linearly Observable Trends

The most important algorithm in ISA (its own description). Finds
A_r ∈ R^(2×q), B_r ∈ R^(q×2), C_r ∈ R^(a×2) minimising the joint
reconstruction error of features and performance:

```
min  ||F~ - B_r Z||_F^2 + alpha*||Y - C_r Z||_F^2   s.t.  Z = A_r F~
```

(the `alpha` weight on the performance term is `opts.pilot.alpha`,
default 1.0, a code addition beyond the original paper's unweighted
objective -- raise it to make PILOT emphasise performance trends over
feature trends in the projection).

- A **global optimum exists but is not unique** (proved in Munoz et al.
  2018, Theorem 1/2); the practical solver is BFGS from `opts.pilot.ntries`
  random restarts (paper default 30; **shipped default is 10**), keeping
  the restart with the highest topological preservation (Pearson
  correlation between feature-space and instance-space pairwise
  distances). A low `ntries` on a hard instance can under-explore the
  optimum landscape -- raise it if PILOT's fit looks unstable across
  reruns.
- An **analytic** fallback exists (`opts.pilot.analytic = true`) via the
  top-2 (or top-3, in 3D) eigenvectors of `[F~; Y][F~; Y]'`, valid only
  when `F~F~'` is invertible; falls back to numerical automatically if
  the feature matrix is rank-deficient. Only applies when
  `opts.pilot.method = 'standard'` (the default).
- PCA is a *provably suboptimal* solution to this same objective
  (Munoz et al. 2018) -- if someone asks "why not just use PCA," this is
  the citable answer, not a stylistic preference.
- **`opts.pilot.method = 'pls'`** (code addition, not in either paper):
  uses Partial Least Squares (`plsregress`) instead of BFGS/analytic --
  maximises covariance between the projection and the performance
  matrix, and does not require the feature matrix to be full column
  rank. Works at 2D or 3D like the standard method.
- **3D extension (ISA3, Simpson et al. 2025)**: `opts.pilot.dims = 3`
  (replaces the legacy boolean `opts.pilot.ISA3D`, still accepted as an
  alias) redefines Z ∈ R^(n×3) and solves the same problem with
  3-column A_r/B_r/C_r; `PILOTviewpoint.m` then solves a second, small
  optimisation for an optimal *viewing rotation* (2 columns) that
  flattens the 3D space for a specific feature/algorithm subset, with an
  orthogonality penalty (default scaling lambda=0.2, low-sensitivity
  below 0.5), storing the result in `model.pilot.viewpoint`.
  `opts.pilot.viewGroups` (a cell array of algorithm index vectors)
  requests one viewpoint per group instead of a single global one. 3D
  retains more information and separates good/bad instances better
  (validated on 56 anomaly-detection variants; F-score improved for 12 of
  16 checked cases moving 2D TRACE3 -> 3D TRACE3), but costs a
  manual/optimised viewpoint choice per subset -- it is not a free
  upgrade for every figure, only where crowding in 2D is the actual
  problem.

### CLOISTER -- Correlated Limits of the Instance Space's Theoretical Regions

Projects the *theoretical* boundary of all feasible instances (not just
the observed ones) by combinatorially enumerating feature-bound vertices
(each feature at its min or max), pruning combinations that violate
observed feature correlations (a vertex cannot pair f_i at its upper
bound with f_j at its lower bound if rho_ij is strongly positive, and the
symmetric rule for strongly negative rho), then taking the convex hull of
the surviving vertices' projections.

- Default correlation-collinearity threshold and significance:
  `opts.cloister.corrThreshold = 0.7`, `opts.cloister.pval = 0.05`
  (matches the paper's epsilon/p notation; the field was renamed from an
  earlier `opts.cloister.cthres`, still accepted as a legacy alias).
- **Code guard beyond the paper**: enumeration is 2^(number of features),
  so the code hard-caps at `opts.cloister.maxFeatures` (default 20) and
  falls back to a plain convex hull of the *observed* projected instances
  if exceeded, with a warning -- at that point CLOISTER is no longer
  estimating a theoretical boundary, it is just re-describing the
  empirical one. Do not interpret a >20-feature CLOISTER boundary as
  theoretical without checking which branch ran.
- The boundary is the tool for assessing benchmark diversity: if real
  instances hug a small region far from the boundary's interior, the
  benchmark suite is not representative of the theoretically possible
  problem space.

### PYTHIA -- Automated Algorithm Selection

Trains one binary good/bad classifier per algorithm on the (z-scored)
2D/3D instance coordinates, predicting PRELIM's binary good/bad label.

- **Code has moved well past the papers' SVM-only description.**
  `opts.pythia.classifier` selects from a registry of MATLAB-native
  classifiers resolved via `ISAgetClassifierFcn`: `'knn'` (`fitcknn`,
  **the default**), `'svm'` (`fitcsvm`), `'tree'` (`fitctree`), `'nb'`
  (Naive Bayes, `fitcnb`), `'linear'` (`fitclinear`), or `'ensemble'`
  (`fitcensemble`). All algorithms in the portfolio share the same
  classifier type. LIBSVM is deprecated for new runs -- `ISAmigrateModel`
  can retrain a legacy LIBSVM-backed model onto the current registry.
- Hyperparameter self-tuning is `opts.pythia.tuning`: `'sobol'` (default;
  scrambled Sobol quasi-random sequence, `opts.pythia.nTuningIter`
  evaluations, k-fold CV via `opts.pythia.kFold`), `'bayes'` (MATLAB
  `bayesopt`, Gaussian-process surrogate, same budget/CV), or `'none'`
  (use pre-supplied `opts.pythia.params` directly, skipping tuning).
  Both `'knn'`/`'svm'`/`'ensemble'` tune two hyperparameters; `'tree'`/
  `'nb'`/`'linear'` tune one.
- Ties among algorithms predicted "good" are broken by higher CV
  precision (`out.Yhat .* trainPrecision`, arg-maxed per instance,
  `computeSelection` in `PYTHIA.m`); if no classifier predicts good for
  an instance, `selection0` (the "confident" recommendation column) is 0
  and `selection1` (the "always give an answer" column) falls back to
  whichever algorithm is most often labelled good across the *training*
  set (`argmax(mean(Ybin))` -- the most frequently-good algorithm, not
  necessarily the one with the single highest continuous performance
  value). Read `selection0` when you want to know which instances PYTHIA
  is actually confident about.
- `opts.pythia.skip = true` bypasses classifier training entirely; TRACE
  then falls back to the true labels directly (with a warning) instead
  of PYTHIA's predictions.

**Precision matters more than accuracy here**: a high-precision,
lower-recall classifier is exactly what you want, because the claim being
made is "when I recommend this algorithm, trust it," not "I never miss a
good instance." Read PYTHIA's precision column before its accuracy
column.

### TRACE / TRACE3 -- Algorithm Footprints

A footprint is characterised by three numbers, all normalised against
the convex-hull baseline of the full instance space: **area** (percentage
of the space -- a robustness proxy), **density** (instances per unit
area relative to the space's overall density -- strength of evidence),
and **purity** (fraction of enclosed instances that are actually good --
absence of contradicting evidence).

- **Legacy TRACE** (Munoz & Smith-Miles 2017; `opts.trace.method =
  'legacy'`, 2D only): DBSCAN to find dense clusters of good instances (k
  and epsilon auto-set from instance count via a density-based
  heuristic), alpha-shapes around each cluster, then pairwise overlap
  resolution between algorithms' best-footprints by removing the
  lower-purity side of any overlap (ties kept, as evidence is
  insufficient to prefer either). `opts.trace.contra` (default true for
  this method only) toggles this contradiction-removal step.
- **TRACE3** (Simpson et al. 2025; **the default method**,
  `opts.trace.method = 'trace3'`, works natively in both 2D and 3D):
  builds an alpha-shape only from instances where PYTHIA's predicted
  label and the true label agree (`Zu = {z_i : yhat_i=1 AND ybin_i=1}`),
  then iteratively shrinks the alpha-radius until a purity threshold
  (`opts.trace.PI`, default 0.6) is met or the shape collapses, trimming
  a fixed fraction of residual area each iteration to discard
  outlier-driven offshoots. Candidate footprints below
  `opts.trace.minInstances` (default 4) or `opts.trace.minAreaFrac`
  (default 0.01, as a fraction of the space) are discarded.
- **The PYTHIA -> TRACE3 coupling is unconditional, not configurable.**
  Earlier code versions exposed `opts.trace.nn`/`opts.trace.prior` fields
  suggesting TRACE3 could train its own from-scratch k-NN classifier
  (matching the ISA3 paper's literal description) as a fallback when
  PYTHIA predictions were unavailable -- **those fields have since been
  removed from `options.json` and `TRACE.m` entirely** (there is no
  `opts.trace.nn`/`.prior` in the current toolkit); TRACE3 always reuses
  PYTHIA's `Yhat` when available, and only falls back to the true labels
  (`Ybin`) directly if PYTHIA was skipped (`opts.pythia.skip = true`).
  If you need the ISA3 paper's literal from-scratch KNN-on-features
  fallback behaviour, it does not currently exist in this toolkit --
  verify against the exact commit in use before assuming otherwise.
- TRACE3 vs. legacy, per the ISA3 paper's own validation: TRACE3 produces
  a footprint (possibly small) in almost every case where legacy DBSCAN
  produces none; TRACE3 footprints have normalised density > 1 far more
  often (denser-than-average regions actually found); legacy sometimes
  reports a higher raw purity, but by covering a much smaller, weak-
  evidence area (excluding real bad instances rather than finding
  where good instances concentrate). Treat legacy TRACE's "no footprint"
  outcome for an algorithm as ambiguous, not as proof of a weak
  algorithm -- this exact failure mode (LOF) contradicted a prior
  published finding and was traced to the footprint algorithm, not to
  LOF actually being weak.
- The purity threshold and footprint-size guards are described in the
  ISA3 paper as low-sensitivity in a wide range (purity threshold 0.6
  default) but each has a documented failure mode at the extremes: too
  permissive admits low-density noise as a "footprint"; too strict
  collapses genuinely good but sparse regions to nothing. Raising the
  purity threshold has little effect unless paired with tighter
  supporting evidence (a higher-precision PYTHIA classifier), not
  adjusted independently.
- **Eval-mode footprints (`explore()`) are re-scored, not rebuilt.**
  Evaluating a trained model on new instances re-tests the *existing*
  footprint polygons against the new points; it does not re-run TRACE3's
  alpha-shape construction. An algorithm present in the new data but
  absent from the trained model gets a placeholder empty footprint, with
  a "not enough instances to calculate a footprint" message -- this is
  expected, not a bug, for a genuinely new algorithm evaluated against an
  existing space.

---

## 4. Reading the results (Step 5 of the methodology)

Work through these in order; each answers a different adequacy question.

1. **Distribution of instance sources**: are real-world instances
   distinct from synthetic ones? Do randomly generated instances cluster
   near the space's centre (typical for naive random generators) while
   structured/real instances occupy distinctive regions? Are there
   visible holes -- regions inside the theoretical boundary with no
   instances at all?
2. **Algorithm-selection accuracy** (PYTHIA precision/recall per
   algorithm): if precision is high and the predicted-good regions
   visually match the empirical good-performance scatter, the selected
   features are adequate. If not, the feature set -- not necessarily the
   algorithms -- is the thing to revisit.
3. **Footprints**: which algorithms have unique, dense, pure regions?
   Overlap and shared strength are real, informative findings, not
   failures of the method -- two strong, broadly competitive algorithms
   with heavily overlapping footprints is itself the finding.
4. **Feature distributions across the space**: for each axis-defining
   feature, does its gradient across the space explain why the easy/hard
   or algorithm-A/algorithm-B regions fall where they do? This is where
   the "why," not just the "where," comes from -- an instance space with
   good footprints but uninterpretable feature gradients has captured
   discriminative power without captured mechanism.
5. **Insights**: synthesise 1-4 into statements a domain expert could
   falsify (e.g. "algorithm X wins when slack is high and teacher-
   conflict degree is low"), not just "algorithm X has R^2 = 0.8."

---

## 5. When to augment vs. stop (Step 6)

Augment the meta-data when: instances are not diverse/dense enough to
fill the theoretical boundary; footprints are near-identical across a
supposedly diverse algorithm portfolio (may indicate the algorithms are
mechanistically similar, or that the instances fail to expose real
differences -- distinguish these two explanations before concluding
either); or PYTHIA's accuracy is poor for specific algorithms (usually a
feature-inadequacy problem, not an algorithm problem). New instances,
new algorithms, and new features are three independent levers -- identify
which one the Step-5 analysis actually implicates before spending effort
generating more of the wrong thing.

Stop when the boundary's interior is adequately filled and an added
iteration selects the same features. Do not treat a single ISA pass as
final by default; say so explicitly if scope or time genuinely limits the
work to one iteration.

---

## 6. Pitfalls (checked against the papers and the live code)

- **The "good performance" threshold is a hidden researcher degree of
  freedom.** State it, justify it, and consider reporting sensitivity to
  it, exactly as you would a significance level.
- **Small feature counts silently change the pipeline.** SIFTED skips
  its clustering+GA stage below 3 or below K features; CLOISTER silently
  degrades to an empirical convex hull above `opts.cloister.maxFeatures`
  (default 20) features. Know which branch executed before interpreting
  the result as the textbook algorithm.
- **PILOT's global optimum is not unique.** Different valid solutions can
  give visually different (but equally optimal by the objective) axes;
  the topological-preservation tiebreak is what actually selects the
  reported figure, and it depends on `opts.pilot.ntries` random restarts
  (shipped default 10, not the papers' 30) -- a low `ntries` on a hard
  instance can under-explore the optimum landscape.
- **A missing footprint is not evidence of a weak algorithm** -- it may
  be a footprint-algorithm artefact (documented LOF/legacy-TRACE case),
  or simply an algorithm not present in the trained model being evaluated
  in `explore()` (placeholder empty footprint, not a real "no evidence"
  result). Cross-check against a second footprint method (legacy vs
  TRACE3, or 2D vs 3D) before concluding an algorithm is uncompetitive.
- **CSV formatting failures are the most common practical error**:
  `NA` instead of `NaN`, Excel error codes (`#REF!`, `#DIV/0!`), or empty
  rows all corrupt PRELIM silently or loudly depending on downstream
  step. Validate the raw CSV before debugging the pipeline.
- **Ground-truth features are legitimate in ISA** (unlike in meta-
  learning/algorithm selection proper) precisely because ISA's goal is
  explanation, not deployable prediction -- do not reflexively flag their
  use as a leakage bug; check which goal the analysis is actually serving
  first.
- **A published ISA figure and a rerun of the current toolkit can
  legitimately differ** even on the same nominal method name (SIFTED and
  PYTHIA are the clearest examples here -- both have moved well past what
  their originating papers describe). When reproducing a specific
  paper's numbers, pin the toolkit commit, not just the method name.
- **Option field names have already changed across recent versions**
  (e.g. `opts.cloister.cthres` -> `.corrThreshold`, `opts.pythia.cvfolds`
  -> `.kFold`, `opts.parallel.flag`/`.ncores` -> `opts.general.parallel`/
  `.ncores`, and several fields removed outright -- `opts.trace.nn`/
  `.prior`, `opts.pythia.uselibsvm`, `opts.trace.usesim`,
  `opts.sifted.NTREES`). A saved `options.json` or `model.mat` from an
  older run may use old names; `ISAvalidateOpts`/`ISAdefaults` accept
  common legacy aliases, and `ISAmigrateModel` brings an old `model.mat`
  forward, but do not assume a config snippet from an older tutorial or
  paper still matches field-for-field.

---

See `references/matlab-toolkit.md` for the `InstanceSpace` MATLAB
repository's file format, `options.json` schema with verified current
defaults, the `InstanceSpace` class API, and pipeline invocation.
