# What's New

Changes in recent versions of the toolbox

The complete list of changes, including every bug fix, is in `RELEASE_NOTES.md` in the repository root.

## v0.9.2 (in development)

### Bug Fixes

- `PYTHIA` evaluation scores each algorithm only on the new instances with observed performance, and reports `NaN` for a trained algorithm that the new data does not cover. Before, such an algorithm was scored against labels that were never measured.
- The Oracle row of the PYTHIA summary reports the fraction of instances on which any algorithm is good. Before, it was always 1, which is wrong when `opts.perf.AbsPerf` is `true`.
- `CLOISTER` computes a 3D boundary for a 3D instance space; `obj.plot('boundary')` and `scriptpng` draw it. Before, the boundary ignored the third coordinate.
- Footprint boundaries written by `scriptcsv` include the boundaries of holes inside a footprint region.

### Documentation

- Reference documentation in the MATLAB Help browser and on the web, with examples for every function.

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
