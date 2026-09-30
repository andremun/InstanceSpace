**InstanceSpace code review — 30 September 2026**

The highest-priority findings are stale SIFTED fitness values, inconsistent portfolio indices after pruning, test-label leakage in selection, and inconsistent training/evaluation transformations. These can change scientific results without an obvious failure. The ordinary reference-data path passed a focused smoke test.

Reviewed checkout: `7fe4b91`. The checkout advanced from `98a01ac` during the review; the intervening changes were inspected and the newly available `doc/` reference pages were included. No production source files were edited. This directory contains the report and diagnostic reproduction artifacts only.

**Scope and evidence**

Reviewed the class API, ingestion and preprocessing, all core pipeline stages, classifier dispatch and migration, output helpers, wrappers, test architecture, and documentation deployment. Context came from `doc/src/`, README, the refactor plan v1.7, Smith-Miles and Muñoz (2023), and Simpson et al. (2025) under `docs/`. The refactor plan is historical: deliberate changes such as PYTHIA/TRACE coupling and the current SIFTED implementation were not treated as defects merely because they differ from an older paper.

“Reproduced” means a focused execution in MATLAB R2026a Update 5. “Static” means verified from the code path, without a dedicated execution. The complete R2025a integration suite, parallel-worker behavior, legacy LIBSVM runtime, and large-data performance benchmarks were not run. No timing speedups are claimed.

P1 means address before relying on affected results. P2 means a concrete correctness, robustness, or supported-workflow defect that should be fixed next.

**Priority findings**

1. **[P1, reproduced] SIFTED's persistent cache survives between calls.** [SIFTED.m:225](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/SIFTED.m:225), [SIFTED.m:250](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/SIFTED.m:250).

   `costfcn` caches fitness by feature-selection bitmask alone. `clear costfcn` inside the sibling `clearCache` function does not reset that local function's persistent state in the tested MATLAB runtime. The first six-feature call evaluated PILOT for six candidates; a second call with changed labels evaluated PILOT zero times. The worker cleanup does not solve the serial-client cache. Repeated experiments can therefore select features using scores from a different dataset or labeling rule.

   **Fix:** create a fresh map per SIFTED invocation and pass it through a nested fitness closure, or implement an explicit reset operation in `costfcn`. Worker caches also need an invocation identity or explicit reset. Test two different datasets/label vectors with identical feature counts against independent clean-session runs.

2. **[P1, reproduced] Algorithm pruning leaves `P` in the original index space.** [InstanceSpace.m:552](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/InstanceSpace.m:552).

   `runPrelim` removes columns from `Yraw`, `Y`, `Ybin`, labels and normalization parameters, but retains `P` calculated before removal. A three-algorithm input with the first algorithm pruned retained two labels while every `P` remained `3`. TRACE's `P==i` tests then lose or misattribute best-performance footprints; portfolio CSVs contain invalid indices. `beta` also retains the original portfolio denominator, which can disagree with the retained portfolio and evaluation-time recomputation.

   **Fix:** define one retained-portfolio mapping and apply it to every algorithm-indexed field. Recompute best labels and easy-instance labels against the intended portfolio. In absolute-threshold mode an algorithm with no good instances can still be best on a hard instance, so merely subtracting deleted indices is insufficient: decide whether pruning itself should be retained, or recompute winners from retained raw performance. Add a test deleting the first and middle algorithms.

3. **[P1, reproduced] Evaluation recommendations use the test labels to choose a fallback.** [PYTHIA.m:452](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/PYTHIA.m:452), [PYTHIA.m:778](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/PYTHIA.m:778).

   `computeSelection` chooses its default from `max(mean(Ybin))`. In evaluation mode, these are the test-set ground-truth labels. Keeping the same trained classifiers and coordinates while swapping test labels changed `selection1` from algorithm 1 to algorithm 2. A newly introduced, untrained algorithm can also win this fallback. The fallback for old models without stored precision additionally derives voting precision from test labels.

   **Fix:** persist the training-time default algorithm and selection weights, and use them unchanged for inference. Use training data to recover these fields during migration where possible; otherwise use a declared deterministic fallback or report them unavailable. Test that changing test outcomes never changes recommendations.

4. **[P1, reproduced] PLS trains and evaluates in different coordinate systems.** [PILOT.m:149](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/PILOT.m:149), [InstanceSpace.m:911](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/InstanceSpace.m:911).

   Training stores `out.Z = XS` from mean-centered PLS, but evaluation uses `X*A'` without subtracting the fitted training mean. This matters when preprocessing is off and also when instance subsetting leaves a nonzero training mean. A direct same-data reproduction had a maximum coordinate discrepancy of about 8.72. Classifiers and footprints trained in one space are applied to shifted points.

   **Fix:** store the PLS feature mean or an affine projection offset. Apply the same transform during training, exploration, and CLOISTER boundary projection. Test that evaluating the training matrix reproduces the saved training coordinates with normalization both enabled and disabled and with row subsetting enabled.

5. **[P1, reproduced] Constant columns become non-finite during evaluation.** [PRELIM.m:192](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/PRELIM.m:192), [PRELIM.m:212](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/PRELIM.m:212).

   Training uses `zscore`, which returned finite zero-normalized values and a stored standard deviation of zero for a constant feature. Evaluation divides directly by that stored zero. On identical data the training matrix was finite and the evaluation matrix was not. The performance-column path has the same direct division. PYTHIA also directly divides by stored coordinate scales at [PYTHIA.m:365](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/PYTHIA.m:365).

   **Fix:** define a constant-column policy at fit time: drop such features with a saved mask, or store an effective scale of one and a constant-column flag. Apply it consistently at evaluation. Test constant features, constant algorithm performance, and degenerate projected dimensions. Do not replace resulting NaNs after the fact; preserve the transform's contract.

6. **[P1, static] Selector statistics labeled as CV use in-sample predictions.** [PYTHIA.m:312](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/PYTHIA.m:312), [PYTHIA.m:766](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/PYTHIA.m:766), [PYTHIA.m:791](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/PYTHIA.m:791).

   Per-algorithm statistics use `Ysub`, but algorithm selection and the Selector row use `Yhat`, the predictions from the final model on its own training data. Selector precision and recall are then presented under `CV_model_precision` and `CV_model_recall`. An overfitting learner can therefore yield optimistic selector statistics in the same table as actual fold predictions.

   **Fix:** calculate and retain separate out-of-fold and fitted-data selections. Build CV selector metrics from out-of-fold predictions, and label resubstitution metrics explicitly. For unbiased end-to-end predictive assessment, use an outer holdout/CV loop that also fits preprocessing, supervised projection, feature selection, and tuning within each training fold; current PYTHIA CV alone does not assess that whole pipeline.

7. **[P1, static] A staged rebuild can freeze options that did not produce the retained artifacts.** [InstanceSpace.m:230](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/InstanceSpace.m:230), [InstanceSpace.m:844](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/InstanceSpace.m:844).

   When a pipeline is complete, `model.opts = obj.opts` copies every current option, even if only TRACE or CLOISTER was rebuilt. For example, changing `opts.norm.flag` and rebuilding TRACE leaves the original preprocessing/projection in place but records a different normalization policy. Exploration then follows that new policy against the old fitted parameters and classifiers. Constructor-only validation also does not validate later public option edits.

   **Fix:** record the options consumed by each stage, compare them at build time, and invalidate affected stages or reject an inconsistent partial rebuild. Persist fitted preprocessing settings with the preprocessing artifact. Revalidate changed options before execution, and make inherited seed/verbosity behavior explicit after object construction.

**Further correctness and robustness findings**

8. **[P2, reproduced] Partial models cannot be loaded.** [InstanceSpace.m:378](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/InstanceSpace.m:378), [InstanceSpace.m:416](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/InstanceSpace.m:416).

   A fresh `build('stages',{'prelim'})`, `save()`, `InstanceSpace.load(...)` fails with `Unrecognized field name "opts"`. `model.opts` is assigned only after every stage completes, whereas public `save()` accepts a partial model and `load()` explicitly attempts to support partial stages.

   **Fix:** serialize stage provenance and required options for partial models, and restore the completed-stage state without pretending missing stages are complete. Add partial save/load/resume tests after each stage.

9. **[P2, reproduced] Re-running SIFTED cannot recover previously removed features.** [InstanceSpace.m:651](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/InstanceSpace.m:651).

   SIFTED overwrites `model.data.X`, labels, and `featsel.idx`. A later SIFTED run receives this reduced matrix. In the reproduction, four features became two, and disabling SIFTED followed by rebuilding that stage still left two. Lowering a threshold or increasing the cluster count likewise cannot revisit discarded features.

   **Fix:** preserve the pre-SIFTED processed input and original row mapping as an immutable stage artifact. Every SIFTED run should start from it; invalidate later stages normally. Test a strict selection followed by a relaxed selection against a clean rebuild.

10. **[P2, reproduced] TRACE unnecessarily computes a test-set convex hull.** [TRACE.m:95](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/TRACE.m:95).

    `convhull(Z)` executes before the evaluation branch, even though that branch uses the trained space geometry. One-point evaluation fails with `NotEnoughPtsConvhullErrId`; collinear 2D or coplanar 3D batches can fail similarly. Membership queries against an existing footprint do not require a new hull.

    **Fix:** move hull construction into training mode. Add one-instance, duplicate-instance, and lower-dimensional evaluation tests in both 2D and 3D.

11. **[P2, reproduced] Relative performance mishandles zero maxima and negative scores.** [PRELIM.m:100](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/PRELIM.m:100), [PRELIM.m:119](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/PRELIM.m:119).

    For maximization with all-zero scores and a 5% tolerance, both algorithms were marked bad, including the best algorithms. The denominator is replaced with `eps`, but `Yaux` retains zero, producing a relative loss of one. Minimizing `[-10,-1]` marked both algorithms good: division by a negative best score reverses the intended comparison.

    **Fix:** define and document the metric domain. Either reject negative scores in relative mode or use a signed-difference formulation with an explicitly defined positive scale. Handle a zero best score as a separate case and preserve raw `Ybest` for reported oracle performance. Test minimization/maximization, zero ties, negative values, and mixed signs.

12. **[P2, reproduced] Supported sparse missing values reach algorithms that cannot consume them.** [INIT.m:243](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/INIT.m:243), [PILOT.m:120](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/PILOT.m:120).

    INIT retains rows with only some missing features, and PRELIM preserves those NaNs. A single missing feature in an otherwise ordinary analytic PILOT input fails at `rank(X)`. Missing performance values can also contaminate the analytic eigendecomposition. In the numerical path, a feature NaN propagates to an entire projected row and ultimately invalid geometry. The documentation says missing cells are allowed, but there is no complete missing-data policy through projection.

    **Fix:** fit and persist a missing-feature treatment, such as a declared imputation rule, or reject incomplete rows at ingestion with their identifiers and a clear policy. Treat unavailable algorithm outcomes as unavailable labels rather than silently equating them with observed failure when reporting classifier quality. Test sparse missing X and Y separately.

13. **[P2, reproduced] Training-time feature removal is not replayed on test metadata.** [INIT.m:140](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/INIT.m:140), [INIT.m:251](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/INIT.m:251).

    A training feature dropped by `nanThreshold` reduces the fitted transform width. Test ingestion does not apply the same drop and instead rejects a file with the original matching column schema. The four-column reproduction trained on three retained features and failed evaluation with `featureCountMismatch`.

    **Fix:** persist retained input feature names/masks and align test columns to them before applying fitted transforms. Preserve checks for genuinely missing required features. Do not independently decide which columns to drop from the test distribution.

14. **[P2, reproduced] Selector recall counts successful selections as false negatives.** [PYTHIA.m:800](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/PYTHIA.m:800).

    `fb = sum(any(Ybin & ~sel0,2))` counts a row whenever any unselected algorithm is good, even when the selected algorithm is also good. With two algorithms good on every instance and perfect selections, reported Selector recall was 50%.

    **Fix:** define selector recall at instance level, using disjoint success/miss events. For example, a missed opportunity is an instance with at least one good algorithm but no good selected algorithm. Add cases with multiple good algorithms, abstentions, and no available good algorithm.

15. **[P2, reproduced] Cost-sensitive weights differ from the documented objective.** [PYTHIA.m:144](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/PYTHIA.m:144), [README.md:202](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/README.md:202).

    README specifies `abs(Y-Ybest)`; the code uses `abs(Y-nanmean(Y(:)))`. On the reproduced first row, the code assigned weights `[9.5,0.5]` where regret-based weights were `[0,10]`. This changes the learning objective substantially.

    **Fix:** implement the intended per-instance regret weights with a documented zero-weight policy, or explicitly revise the public contract if deviation from the global mean is the intended objective. Add a numeric weight assertion; the current option test only checks that training runs.

16. **[P2, reproduced] The 3D camera azimuth uses the wrong convention.** [PILOTviewpoint.m:153](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/PILOTviewpoint.m:153), [scriptfcn.m:192](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/output/scriptfcn.m:192).

    `cart2sph` measures azimuth from positive X; the stored value is passed directly to MATLAB `view`. The runtime reproduction requested direction `[1,0,0]` through this conversion but obtained camera direction `[0,-1,0]`. The optimized plane is therefore not the plane shown by the output camera in general.

    **Fix:** store a view direction and use MATLAB's vector view API, or explicitly convert between spherical and camera azimuth conventions. Verify rendered camera direction against the cross product for basis vectors and non-axis-aligned views, including aspect-ratio handling.

17. **[P2, reproduced] Equal-length viewpoint groups fail the JSON round trip.** [ISAvalidateOpts.m:282](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/utils/ISAvalidateOpts.m:282), [PILOTviewpoint.m:88](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/PILOTviewpoint.m:88).

    `jsondecode(jsonencode(...))` turns `{[1 2],[3 4]}` into a numeric matrix. Validation requires a cell array and rejects it before PILOTviewpoint's numeric-matrix compatibility code can run. The existing test deliberately uses unequal group lengths, avoiding this failure.

    **Fix:** canonicalize decoded groups before validation, or accept and validate the numeric representation at the public boundary. Test equal-length, unequal-length, singleton, and empty groups.

18. **[P2, reproduced] Valid master seeds overflow PYTHIA's derived seed.** [PYTHIA.m:504](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/PYTHIA.m:504).

    `foldSeed = baseSeed*1e5 + fold*1e3` can exceed MATLAB's supported seed range. Seed 50000 failed with `MATLAB:rng:badSeed`, although it is a valid master seed.

    **Fix:** use bounded deterministic substreams or a specified integer seed derivation, and validate the actual master-seed range. Check boundary seeds and consistency across serial and parallel execution.

19. **[P2, reproduced] SIFTED's analytic fallback lacks numerical options.** [SIFTED.m:233](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/SIFTED.m:233), [PILOT.m:120](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/PILOT.m:120).

    SIFTED calls PILOT with `analytic`, `dims`, and `verbose`, but no `ntries`. If the selected feature matrix is rank-deficient, PILOT switches to its numerical path and accesses `opts.ntries`. The matching standalone invocation failed with `Unrecognized field name "ntries"`.

    **Fix:** populate complete PILOT defaults at the standalone boundary, including defaults for fallback paths. Pass the configured SIFTED seed into the inner projection as well. Test collinear and constant candidate subsets.

20. **[P2, reproduced/static] Option normalization and subsetting dispatch disagree.** [InstanceSpace.m:132](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/InstanceSpace.m:132), [InstanceSpace.m:581](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/InstanceSpace.m:581), [ISAvalidateOpts.m:260](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/utils/ISAvalidateOpts.m:260).

    A scalar-string root directory without a trailing slash becomes a two-element string array under bracket concatenation and fails the `rootdir` property assignment (reproduced). A valid string `selvars.type="Ftr"` passes validation but fails the `ischar` gate and silently disables density filtering (reproduced). Case-insensitive tuning validation accepts `BAYES`/`NONE`, but PYTHIA's case-sensitive dispatch can run Sobol instead (static). Enabling both density and fractional subsetting can reference density outputs that were never computed (static).

    **Fix:** normalize paths with a single char/string policy and `fullfile`, canonicalize enum values once, and validate mutually exclusive modes. When file-index filtering is requested, reject missing files and invalid indices instead of silently switching to the full dataset.

21. **[P2, static] Repeated exports can leave obsolete footprints looking current.** [scriptcsv.m:54](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/output/scriptcsv.m:54), [scriptfcn.m:452](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/output/scriptfcn.m:452).

    Footprint CSVs are only written when the current boundary is nonempty. A rerun yielding no footprint leaves the old file in place. A subsequent 3D run also skips boundary extraction entirely, so old 2D footprint CSVs can survive beside new 3D coordinates. Algorithms removed from a later portfolio leave their earlier outputs behind too.

    **Fix:** write outputs into a fresh run directory and publish a manifest, or remove only files owned by the previous run before replacing them. Represent empty/unsupported geometry explicitly. Add a two-run test in the same output directory.

22. **[P2, static; already acknowledged in repository notes] Geometry export remains incomplete.** [CLOISTER.m:76](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/CLOISTER.m:76), [CLOISTER.m:106](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/CLOISTER.m:106), [scriptfcn.m:496](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/output/scriptfcn.m:496).

    CLOISTER computes only an XY convex hull even for 3D projections; `scriptcsv` can still export those selected vertices as three-column bounds. `traceOneRegion` traces only one boundary cycle of an alpha-shape region, omitting holes. Its starting cycle is arbitrary, so the implementation does not itself establish that the retained cycle is the outer boundary. Three-dimensional footprint surface CSV export is absent.

    **Fix:** represent 3D hulls as vertices plus face connectivity, trace every 2D boundary cycle with explicit hole semantics, and test reconstructed geometry against the stored shape. Until supported, identify unavailable geometry in the exported schema rather than supplying a misleading approximation. CLOISTER 3D and hole tracing are already noted as deferred issues #50/#52 locally; live issue status was not checked.

23. **[P2, static] The documentation deployment omits its stylesheet and landing index.** [docs-pages.yml:25](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/.github/workflows/docs-pages.yml:25), [generate.py:55](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/doc/generate.py:55).

    The deployment uploads only `doc/html`, but all generated pages reference `../style.css`, located in `doc/style.css`, outside the artifact. The artifact also contains `Landing.html` rather than `index.html`. The root deployment URL has no landing index and the page stylesheet cannot be served from the uploaded artifact. The workflow uploads committed HTML without regenerating it.

    **Fix:** generate a self-contained site directory containing the HTML, stylesheet, and an index page, then publish that directory. Validate local links and artifact contents in CI. Correct nearby documentation contracts too: `pythia.flag` is advertised as enabling/disabling the stage but is not consumed, `betaThreshold` is a fraction of algorithms rather than instances, and several value-class walkthrough calls omit reassignment to `obj`.

**Behavior requiring an explicit design decision**

- **TRACE can return a footprint below `PI`.** [TRACE.m:278](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/TRACE.m:278). A reproduced shape returned purity 0.5 with `PI=0.9`. Spectrum exhaustion is an allowed stopping condition in the historical plan, so this is not automatically an implementation departure. However, the public description calls `PI` a minimum acceptable purity. Decide whether unsuccessful candidates should be empty, or retain them with `accepted=false` and a termination reason. The current “best footprint found” comment is also misleading: no best-so-far candidate is tracked.
- **Failed CV folds can become ordinary reported results.** [PYTHIA.m:518](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/PYTHIA.m:518), [PYTHIA.m:650](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/PYTHIA.m:650). When all Sobol candidates fail, the first candidate is selected anyway. Fixed-parameter tuning does not reject the NaN fold scores at all. If final fitting succeeds, placeholder false predictions feed normal-looking statistics. Fail explicitly or mark the model and affected metrics invalid; handle single-class folds as a defined case.
- **Score columns are not always probabilities.** SVM posterior-fit failure deliberately falls back to raw decision scores, while `Pr0hat` remains documented as probability. One-class fold scores also should be interpreted through `ClassNames`, not blindly through column 1. Expose score type/calibration status and map labels explicitly.

**Performance improvements, in recommended order**

| Location | Current cost | Proposed improvement and trade-off |
|---|---|---|
| [PILOT.m:108](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/PILOT.m:108) | Always computes all pairwise feature distances, including analytic/PLS/precalculated paths that never read them. Profiling confirmed the analytic path enters `pdist`. A single double distance vector for 20,000 rows is about 1.60 GB decimal (1.49 GiB). SIFTED repeats this for candidate subsets. | Move distance construction into the numerical restart-ranking branch. This removes unnecessary work without changing results. For numerical runs, consider blockwise exact correlation accumulation or a fixed sampled pair set; sampling changes the ranking objective and must be validated. |
| [CLOISTER.m:83](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/CLOISTER.m:83) | Enumerates `2^p` corners and performs nested feature-pair checks. At the permitted default limit of 20 features, `idx` and `Xedge` alone occupy about 320 MiB, excluding temporaries; worst-case pair checks approach 199 million. | Precompute only significant correlation constraints, process corners in bounded batches, and avoid retaining both full arrays. Use an explicit memory budget. Any replacement of the empirical bound by an observed-data hull changes scientific meaning and should remain an explicit fallback. |
| [PYTHIA.m:690](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/PYTHIA.m:690) | Fits an SVM posterior for every candidate in every outer CV fold, though candidate selection uses label error. Posterior fitting introduces its own additional fitting workload. | Tune using classification labels/scores, then calibrate the selected model and any outputs that actually require calibrated probabilities. Verify that the selected labels and calibration contract remain equivalent. |
| [TRACE.m:325](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/TRACE.m:325) | Queries shape membership for all points, then again for the good subset at each shrink step. | Compute `inside = inShape(as,Z)` once and obtain good count with `sum(inside & Ybin)`. Consider evaluating actual critical alpha values rather than a 100-point linear grid, but treat that as an algorithm change requiring geometry comparison. |
| [PILOT.m:224](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/PILOT.m:224), [PILOTviewpoint.m:135](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/PILOTviewpoint.m:135) | Numerical finite differences repeatedly evaluate large residual matrices; viewpoint search allows 30,000 iterations and a `1e-20` function tolerance. | Supply verified objective gradients, cache invariant products where appropriate, and select tolerances using measured reconstruction/view stability. Record convergence status rather than assuming all restarts converged. |
| [FILTER.m:108](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/FILTER.m:108) | `rangesearch(X,X,r)` materializes every neighborhood. It avoids a dense distance matrix for sparse neighborhoods, but duplicates/dense data still give quadratic stored neighbor indices. | Query kept points or batches against a reusable spatial searcher. Preserve the current greedy row order to retain identical selected instances. |
| [SIFTED.m:156](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/SIFTED.m:156) | Computes silhouette diagnostics for every K through the total feature count even though only `opts.K` is used for selection. | Make the advisory sweep optional or bound it around the requested K. Separate diagnostics from required selection work. Fix cache correctness before optimizing cache hit rate. |
| [InstanceSpace.m:507](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/InstanceSpace.m:507) | Pool startup/teardown occurs per staged build when the build owns the pool, and errors bypass normal cleanup. Core stages may use an existing pool even when `general.parallel=false`. | Define explicit pool ownership, use `onCleanup`, and pass a resolved execution policy into stages. Avoid deleting unrelated user pools. Measure pool overhead separately from algorithm time. |

The analytic projection also forms `X'*X` before solving at [PILOT.m:188](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/core/PILOT.m:188). A QR/SVD-based least-squares solve avoids squaring the condition number for nearly collinear features. This is primarily a numerical-stability improvement; test reconstruction and projection equivalence before replacing it.

**Validation and regression coverage**

The diagnostic drivers are [isa_review_checks.m](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/review/2026-09-30/isa_review_checks.m), [isa_review_extra.m](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/review/2026-09-30/isa_review_extra.m), and [isa_review_smoke.m](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/review/2026-09-30/isa_review_smoke.m). They use synthetic data or copies of the reference CSVs in temporary directories. They deliberately catch and print expected failures; they are diagnostic reproductions, not a passing/failing `matlab.unittest` suite. They add the reviewed checkout to the MATLAB path. Raw results are in [reproductions.log](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/review/2026-09-30/reproductions.log) and [smoke.log](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/review/2026-09-30/smoke.log).

The smoke configuration used analytic 2D PILOT, SIFTED disabled, KNN with two Sobol candidates, absolute error threshold 0.2, and image/CSV output disabled. Build, save, load, and exploration completed: 212 training rows, 235 test rows, ten classifiers, finite train/test coordinates. This does not establish correctness of numerical/3D/parallel or export paths.

The current option-coverage suite calls build/explore and succeeds if no exception is thrown; it does not assert numerical invariants ([PipelineOptionsTest.m:90](/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace/tests/PipelineOptionsTest.m:90)). That explains why label leakage, wrong weights, shifted projections, and incorrect summary values can coexist with broad option coverage. Targeted tests should establish:

- Training/evaluation transform equality on identical data, including constant columns, PLS, and row/feature subsetting.
- Selector invariance to test labels, and correct instance-level metrics with multiple good algorithms.
- Feature-selection independence across calls and equivalence of staged rebuilds to fresh builds.
- Valid indices and aligned fields after pruning; explicit unknown-label behavior.
- Save/load/resume of partial models and consistent per-stage option snapshots.
- One-row/degenerate evaluation, sparse missing values, JSON round trips, and seed boundaries.
- Geometry round trips for holes/3D, and removal or invalidation of stale export files.

Fix the P1 state, transformation, and selection problems first. Then add the focused regression cases before performance changes; otherwise faster execution can preserve or amplify incorrect results.
