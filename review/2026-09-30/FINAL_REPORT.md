# Code review resolution report

The review findings were rechecked against `master` at `08044f5`. Changes are on the local branch `codex/fix-sifted-cache-reset`. No push, pull request, or merge was performed. The original `REPORT.md` is retained as a historical record of the earlier checkout, not a description of the final branch.

## Resolution of the 23 findings

| Finding | Resolution | Main commit(s) |
|---|---|---|
| 1. SIFTED cache survives between calls | Explicitly reset client and worker caches before each run and wait for worker resets. Serial and two-worker regression cases compare repeated runs with fresh runs. | `fdbeaa6` |
| 2. Stale portfolio indices after pruning | Recompute preprocessing from retained raw outcomes. Winners, best performance, binary labels, beta, and transforms describe the same retained portfolio. Both passes seed explicitly from `general.seed`. | `2230225`, `c36fa7b` |
| 3. Test-label leakage in recommendations | Persist training fallback and precision weights. Migration recovers them from training data where possible. Legacy models without them use deterministic, training-only defaults. | `d351f99` |
| 4. PLS coordinate mismatch | Store the fitted feature mean and use it for exploration and CLOISTER projections. Recover means for older PLS models when training features are available. | `dcd11e5` |
| 5. Constant-column evaluation failures | Use a unit divisor for constant feature, performance, and classifier-coordinate columns, including old saved zero scales. | `0d5dc64` |
| 6. Mixed training/CV/test summary labels | Retain fitted selections for training plots and footprints. Build the training summary from out-of-fold selections and label exploration metrics as test metrics. | `3ed28bf`, `33f67c2` |
| 7. Staged option provenance | Store options consumed by each stage. Revalidate public option edits and reject partial rebuilds that would retain artifacts fitted with different options. | `5ce3859` |
| 8. Partial save/load | Persist options and stage provenance after each stage. Partial models can be saved, loaded, and resumed without claiming later stages completed. | `5ce3859` |
| 9. Irrecoverable SIFTED feature removal | Preserve the pre-selection input and restart every SIFTED run from it. Legacy models without it must rebuild preprocessing. | `402dde9` |
| 10. TRACE test hull failure | Construct hulls only during training. Evaluation accepts single, repeated, and lower-dimensional test batches. | `18d288a` |
| 11. Zero and signed relative scores | Reject negative raw performance. Restore the original ratio-minus-one / one-minus-ratio expressions and documented zero substitution. Preserve raw best performance; use the same transformed values for labels so zero ties remain good. | `817483c`, corrected by `2c2c36d` |
| 12. Missing-data propagation | After existing ingestion filtering, reject incomplete training features/outcomes with instance identifiers. Reject non-finite projection inputs at the standalone boundary. Missing test outcomes remain unobserved, including rows with no outcomes. | `3ad1174`, `d4e1244` |
| 13. Training feature exclusions not replayed | Save excluded input feature names and remove those columns during test ingestion before checking the retained schema. | `3ad1174` |
| 14. Selector recall double counting | Use disjoint instance-level successes and missed opportunities. A successful choice is not a miss when another algorithm is also good. | `1b45e3b` |
| 15. Incorrect cost-sensitive weights | Use per-instance regret, with the smallest positive regret replacing zeros and uniform weights when all regrets are zero. | `edd7ee4` |
| 16. Camera convention mismatch | Convert spherical azimuth to MATLAB camera azimuth, use equal data-axis scales, and apply the same convention when recalling saved views. | `6945f95`, `61096b1` |
| 17. JSON viewpoint groups | Canonicalise rectangular numeric JSON arrays to one group per row before validation. Test equal, unequal, singleton, and empty group lists. | `6945f95` |
| 18. Seed overflow | Master already bounded fold seeds. Bound algorithm-derived seeds too and validate the supported master/stage seed range. | `6945f95` |
| 19. Incomplete numerical fallback options | Fill standalone PILOT defaults and pass the SIFTED seed into candidate projection. Rank-deficient candidate projection now reaches a complete numerical configuration. | `971ccd1` |
| 20. Inconsistent paths/options/subsetting | Normalise scalar-string directories and enum case. Enforce exclusive subsetting modes and reject missing index files or invalid row indices. | `19a8167` |
| 21. Stale exports | Replace toolkit-owned geometry on each export, including removed algorithms and empty footprints. Plot export also removes obsolete toolkit plot and FIG names. | `c2ac52d`, `ee9a58a` |
| 22. Incomplete geometry | Master already fixed 3D CLOISTER hulls and tracing of all 2D boundary cycles. Add vertex/triangle CSV export for 3D footprints and bounds, explicit empty geometry, and a geometry manifest. | `c2ac52d` |
| 23. Documentation deployment | Already resolved on master: self-contained assets, landing index, generated-source checks, and link validation. Updated behavioral contracts and regenerated the reference site. | Master; `d20aab1` |

Finding 10 moves the hull calculation into training; it does not remove normalization. Exploration retains the training-space measure and density. Footprint area is divided by training-space area, while test-point density is divided by training-space density. Evaluation density therefore depends on test-batch size. The removed evaluation hull calculation was unused and could fail for small or degenerate batches.

Finding 6 is a reporting issue under the intended build/training and explore/testing workflow. These classifier CV summaries remain conditional on the fitted preprocessing, feature selection, projection, tuning, and selection weights. They are not an unbiased outer-CV assessment of the complete ISA pipeline.

## Additional decisions from the review

- TRACE3 keeps its existing alpha-search stopping rule. A candidate below the requested purity is retained with `accepted=false` and `terminationReason='spectrumExhausted'`. Empty footprints carry an explicit reason. Exploration preserves these training acceptance fields while updating test metrics.
- PYTHIA rejects an invalid selected CV result or an entirely failed candidate search. A single-class training fold predicts its observed class. Failed folds cannot silently become ordinary reported statistics.
- Class scores are mapped through class names. `scoreType`, `scoreTypeCV`, and `Pr0subIsProbability` distinguish probabilities, raw scores, and mixed CV outputs. Missing outcomes are excluded from scoring rather than treated as observed failures.

## Performance work and trade-offs

Implemented changes preserve the existing scientific objectives:

- PILOT computes pairwise feature distances only when numerical restart ranking needs them. Analytic projection solves least squares directly instead of using normal equations.
- CLOISTER enumerates corners in batches of 4096 and retains hull vertices between batches. This avoids allocating all feature corners and bit indices at once. The hull remains the hull of the permitted corners, not an observed-data approximation.
- FILTER queries a reusable spatial index one retained instance at a time. It preserves greedy row order without storing every neighbourhood.
- TRACE queries membership once per candidate and derives the good count from that mask.
- SIFTED can omit its advisory silhouette sweep with `sifted.diagnostics=false`; the default retains the diagnostic output.
- Builds respect serial execution settings, reuse existing user pools, and close only their own pool on success or failure.

The proposed SVM calibration shortcut was rejected after an explicit check: posterior calibration changed predicted labels in the synthetic comparison in `check_svm_calibration_labels.m`. Fold calibration is retained. Removing it would change tuning decisions, not merely runtime.

Sampled distance objectives, a different alpha schedule, numerical objective gradients, and relaxed viewpoint tolerances were not adopted. They remain optional optimisation work requiring separate numerical and performance benchmarks. No measured speedup is claimed here. Keeping pre-SIFTED data increases model storage; recomputing PRELIM adds work only when algorithms are pruned.

The relative-performance change in `817483c` was revisited at the user's request. For positive best scores its difference-over-best expression was algebraically equivalent to the original ratio expression, but changed floating-point behavior at epsilon and clipped positive denominators smaller than machine epsilon. The final implementation retains the original expressions directly. The user selected nonnegative raw performance with zero allowed; negative inputs now raise an error instead of using the proposed signed extension. The paper's Algorithm 1 supplies the ratio formulas, and the refactor plan discusses the zero-substitution convention. The finite zero surrogate remains unit-sensitive and is documented accordingly.

## Compatibility and result changes

Affected experiments should be rebuilt to obtain corrected feature selection, portfolio labels, projections, selections, weights, and summaries. Missing training values now raise actionable errors instead of reaching algorithms that cannot consume them. Raw performance is required to be nonnegative; the proposed signed-score extension was withdrawn after the mathematical review. The original relative formulas and zero substitution are retained. Empty CSV geometry and 3D face tables extend the export contract, described by `geometry_manifest.json`. User-authored files must not use the documented toolkit output namespaces.

Training/test separation remains intact: exploration applies saved preprocessing, projection, selection policy, classifiers, and footprints. Test outcomes are used for evaluation metrics, never to choose recommendations.

## Verification

Focused regressions passed in MATLAB R2026a Update 5, including serial/parallel cache isolation, all 16 portfolio-pruning combinations, PLS train/explore agreement with preprocessing on/off and row subsetting, constant columns, label-independent selection, CV summaries, missing-data guards, feature recovery, all-stage save/load/resume, degenerate test geometry, performance-domain checks, positive-score threshold compatibility, seed bounds, JSON groups, camera direction, exact geometry exports, bounded-memory hull equivalence, and pool cleanup on errors.

The documented `example` completed successfully. The complete integration suite passed **198/198 tests**, with zero failed or incomplete tests, in four disjoint clean MATLAB sessions: 177 non-option tests and three groups of seven pipeline-option tests. The helper `run_review_verification.m` runs the same `TestSuite.fromFolder` inventory used by `test_integration`; it partitions the suite to run independent cases concurrently.

After the final relative-performance correction, the affected NumericReviewTest, PortfolioPruningTest, StageUnitTest, StateReviewTest, and ReviewFixTest suites were rerun, followed by the final numeric checks: **80/80 test executions passed** (including repeated numeric cases). The current suite contains two additional numeric tests. Two further PRELIM evaluation regressions and their fresh reference-model build completed without reported test failures; a result-array concatenation error occurred afterward in their reporting command. Those two tests are not included in the structured follow-up pass count. Corrected result collection completed successfully. `verification.json` records every counted test name, status, and duration, and identifies the code revisions.

The first combined example/integration attempt was stopped after source changes during the run made its cached MATLAB code inconsistent with the tests. Its partial output is not used as verification evidence. The four clean-session runs above replaced it.

Documentation generation and validation: `python3 doc/generate.py --check` passed; `.github/scripts/check_doc_links.py doc/html` checked 31 pages with zero problems. The MATLAB Help/Pages output was regenerated and committed.

R2025a is the repository's CI target. Local execution used R2026a; no claim is made that remote CI ran. No GitHub state was changed.
