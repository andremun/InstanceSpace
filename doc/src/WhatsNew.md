# What's New

Changes in recent versions of the toolbox

The complete list of changes, including every bug fix, is in `RELEASE_NOTES.md` in the repository root.

## v0.9.2

### New Functionality

- Independent `INIT` stage for metadata ingestion, callbacks and saved input snapshots. PRELIM reruns reuse the snapshot; rerunning INIT reloads metadata and invalidates downstream results.
- Explicit `opts.trace.boundaryTolerance`, a nonnegative distance in projected coordinates, defaults to zero. Footprints include their exact boundaries; the saved tolerance is reused during exploration.
- Serial GNU Octave 2D/3D build, save/load and exploration with KNN and none/Sobol tuning, TRACE3, CSV and PNG. Validation uses Octave 11.3, Statistics 2.0.0 and Datatypes 1.5.0. Qt and a display are required for 3D PNG. Bayesian tuning, parallel execution, other classifiers, legacy TRACE, FIG and web palette output remain MATLAB-only. Octave classifier archives are runtime-specific.
- 3D CLOISTER boundaries and CSV triangle-face exports, including a geometry manifest and explicit empty footprints.

### Correctness and Reporting

- Recompute preprocessing after pruning algorithms that are never good, reset SIFTED caches between runs, retain inputs for feature-selection reruns, and apply fitted PLS centering during exploration.
- Recommendations depend on training policy and new features. Missing test outcomes are excluded from metrics. Oracle probability uses only rows with observed outcomes; selector recall avoids double counting. Classifier summaries distinguish CV from test metrics and probabilities from raw scores.
- Persist stage options and provenance for safe partial build/save/load/resume; reject incompatible partial rebuilds and replay training feature exclusions during exploration.
- Evaluate footprints against trained geometry for small or degenerate test batches. Preserve holes in 2D exports and record footprint acceptance and search termination. Replace obsolete toolkit-owned output files on export.
- Reproducible seeding immediately before SVM posterior calibration, complete numerical projection defaults, constant-column scaling and validation of seed bounds and input data.

### Upgrading and Limitations

Rebuild experiments from the original metadata to obtain corrected results. Loading an existing model does not apply fixes to its fitted artifacts. Back up models and outputs first; rebuild CLOISTER for older 3D models and preprocessing before rerunning SIFTED if the saved feature-selection input is absent. Raw performance must be nonnegative (zero is allowed), and incomplete training data remaining after ingestion filtering is rejected. See [Migrating a Legacy Model](MigratingLegacyModel.html) for pre-v0.9.0 layouts.

Classifier CV summaries are conditional on the fitted preprocessing, feature selection, projection and tuning. They are not an outer-CV estimate of the full pipeline. Numerical PILOT restart selection can amplify floating-point differences across BLAS implementations and thread settings, so a seed alone does not promise identical projections across environments. `opts.pilot.topoWeight` remains reserved and disabled.

### Documentation

Reference documentation is available in the MATLAB Help browser and on the web, with examples for every function. The version header, release notes and citation metadata now identify v0.9.2.

## v0.9.1

Engineering release; results of the pipeline are unchanged.

- `build` and `explore` accept an `'onStage'` callback, called after each stage.
- `obj.plot('boundary')` and `distribution_boundary.png` show the CLOISTER boundary of a 2D space.
- `INIT` reads the metadata for both training and evaluation, and `PRELIM` gained an evaluation mode, so evaluation now breaks ties exactly as training does.
- `opts.general.seed` now also controls `PILOT` and `SIFTED`.
- `InstanceSpace` checks that each stage's inputs exist before running it.
- Continuous integration and a `matlab.unittest` test suite.

## v0.9.0

Complete refactor of the toolbox.

- The `InstanceSpace` class is the main interface; `buildIS` and `exploreIS` remain as wrappers.
- 3D instance spaces (`opts.pilot.dims = 3`) with optimised viewpoints (`PILOTviewpoint`).
- Partial Least Squares projection (`opts.pilot.method = 'pls'`).
- Classifier registry for `PYTHIA` (`opts.pythia.classifier`) with Sobol or Bayesian tuning, replacing LIBSVM.
- TRACE3 footprints for 2D and 3D spaces.
- `ISAmigrateModel`, `ISAvalidateOpts` and `ISAdefaults`.
- **Licence change**: from GPL v3 to the PolyForm Noncommercial License 1.0.0.

## See Also

[Migrating a Legacy Model](MigratingLegacyModel.html) | [Deprecated Functions](Deprecated.html)
