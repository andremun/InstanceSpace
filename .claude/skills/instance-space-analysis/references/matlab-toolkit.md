# MATLAB `InstanceSpace` Toolkit -- Operational Reference

Verified directly against this repository (`andremun/InstanceSpace`,
v0.9.0). Requires **MATLAB R2025a or later**, with the Global
Optimization, Parallel Computing, Optimization, and Statistics and
Machine Learning toolboxes. Re-verify against the live repo if a detail
below matters for a specific decision -- option names and defaults have
already changed across recent versions (see `RELEASE_NOTES.md` at the
repo root for the version-to-version breaking changes).

**Licence**: PolyForm Noncommercial License 1.0.0 (changed from GPL v3 as
of v0.9.0) -- non-commercial research/educational use only without prior
written permission; see `LICENSE` and `CITATION.cff`. Do not assume the
older GPL terms apply to a current checkout.

## Directory layout

```
InstanceSpace/
  InstanceSpace.m   Primary API: the InstanceSpace value class (build/explore/save/load)
  buildIS.m         Thin backward-compatibility wrapper around InstanceSpace.build()
  exploreIS.m       Thin backward-compatibility wrapper around InstanceSpace.explore()
  startup.m         Adds core/, output/, utils/, deprecated/ to the MATLAB path
  Contents.m        MATLAB Central File Exchange version/date metadata
  CITATION.cff      Machine-readable citation metadata (GitHub/Zenodo)
  RELEASE_NOTES.md  Version-to-version changes, including breaking option renames
  example.m         Minimal getting-started example (single run, bundled dataset)
  test_integration.m  Exhaustive option-coverage regression suite
  liveDemoIS.m      Interactive, stage-by-stage walkthrough (open in MATLAB's Live Editor;
                     a %%-sectioned plain .m file, not a .mlx)
  options.json      A *generated* artifact written by example.m/test_integration.m from
                     the opts struct actually used -- not a source file to hand-edit for
                     lasting effect (the next run overwrites it)
  core/             PRELIM, SIFTED, PILOT, PILOTviewpoint, CLOISTER, PYTHIA, TRACE,
                     TRACE_legacy, FILTER
  output/           scriptcsv, scriptpng, scriptweb, scriptfcn, ISArecallView
  utils/            ISAdefaults, ISAvalidateOpts, ISAgetClassifierFcn, ISAmigrateModel,
                     ISAsubsetData
  deprecated/       PYTHIA2, PYTHIAtest, SIFTED2 -- warn-and-forward shims kept for
                     backward compatibility; do not build new code on these
  test/data/        Bundled reference dataset (metadata.csv, metadata_test.csv,
                     workspace_test.mat) plus test_integration.m's per-case output
                     subdirectories (gitignored except the three reference files)
```

`InstanceSpace.m`/`buildIS.m`/`exploreIS.m` add `core/`, `output/`,
`utils/`, and `deprecated/` to the MATLAB path automatically the first
time any of them is used in a session. To call a function from one of
those folders directly (`PILOT`, `SIFTED`, ...) without going through one
of those three entry points first, run `startup.m` (`run('startup.m')`,
or just `startup` with the repo root as the current folder).

## Required input files (in the working/root directory)

- **`metadata.csv`**: one row per instance. Columns:
  - `Instances` (any case variant) -- unique instance identifier
    (string; numeric IDs are auto-converted to strings). Mandatory.
  - `Source` (optional) -- instance class/family, read as a categorical;
    used for the "distribution of sources" plot.
  - `feature_<name>` -- one column per feature (numeric, `double`). Must
    have more than two features. The `feature_` prefix (case-insensitive)
    is how the pipeline identifies feature columns and is stripped before
    use; multi-word names separated by `_`, no spaces.
  - `algo_<name>` -- one column per algorithm's performance (numeric,
    `double`), same prefix convention. Multiple algorithms in one file.
  - Column order does not matter; only the prefixes and the `Instances`/
    `Source` names do.
  - Avoid `NA` (use `NaN`), Excel error strings (`#REF!`, `#DIV/0!`), and
    empty rows -- these make `readtable` infer a text column instead of
    numeric, which crashes the pipeline downstream rather than failing
    clearly at load time. This is the most common real-world breakage.
- **`options.json`** (optional): configuration; omitted fields are filled
  by `ISAdefaults.m`. Every user-supplied field is validated by
  `ISAvalidateOpts.m` before defaults are filled in, raising a clear
  `ISA:ISAvalidateOpts:*` error on the first invalid value.
- **`metadata_test.csv`** (optional, for `explore()`/`exploreIS.m` only):
  same schema as `metadata.csv`, evaluated against an already-built
  model. Feature columns must match the trained model's feature set in
  both count and order (checked explicitly; mismatches error rather than
  silently misaligning). Algorithm columns are reconciled by name against
  the trained model's algorithm list: known algorithms line up
  positionally; algorithms present here but absent from training are
  appended as new columns and evaluated with an empty/placeholder
  footprint and no trained classifier (not an error).

## `options.json` schema (verified current defaults)

```json
{
  "general":  {"seed": 42, "verbose": true, "parallel": false, "ncores": 18},
  "perf":     {"MaxPerf": false, "AbsPerf": false, "epsilon": 0.05, "betaThreshold": 0.55},
  "prelim":   {"iqrMultiplier": 5, "nanThreshold": 0.20},
  "auto":     {"preproc": true},
  "bound":    {"flag": true},
  "norm":     {"flag": true},
  "selvars":  {"smallscaleflag": false, "smallscale": 0.3,
               "fileidxflag": false, "fileidx": "",
               "densityflag": false, "mindistance": 0.1, "type": "Ftr&Good"},
  "sifted":   {"flag": true, "rho": 0.1, "pval": 0.05, "K": 10,
               "MaxIter": 1000, "Replicates": 100},
  "pilot":    {"analytic": false, "ntries": 10, "dims": 2, "method": "standard", "alpha": 1.0},
  "cloister": {"pval": 0.05, "corrThreshold": 0.7, "maxFeatures": 20},
  "pythia":   {"classifier": "knn", "tuning": "sobol", "nTuningIter": 20, "kFold": 5},
  "trace":    {"method": "trace3", "PI": 0.6, "minInstances": 4, "minAreaFrac": 0.01},
  "outputs":  {"csv": true, "png": true, "web": false}
}
```

`opts.selvars.feats` / `opts.selvars.algos` (not shown above -- absent by
default) restrict the analysis to a hand-picked cell array of feature/
algorithm names, **with** their `feature_`/`algo_` prefix.

Field notes beyond their names -- the full annotated reference is
README.md's "Options" section (kept in pipeline execution order); the
highlights most likely to surprise someone reading only the JSON:

- `perf.MaxPerf`/`AbsPerf` select which of the four PRELIM good-
  performance definitions applies (maximise vs minimise; absolute
  threshold vs relative-to-best margin) -- see SKILL.md Section 3.
  `perf.betaThreshold` sets the fraction of algorithms that must be
  "good" for an instance to count as beta-easy (TRACE's hard footprint).
- `selvars.smallscaleflag`/`fileidxflag`/`densityflag` are three
  mutually-relevant (not combinable in one call) ways to subset instances
  before the pipeline runs: random fraction, an external index file, or
  `FILTER`'s density/representativeness filter (`type` trades off
  feature-space vs algorithm-performance vs good/bad-label similarity
  when deciding which near-duplicate instances to drop: `'Ftr'`,
  `'Ftr&AP'`, `'Ftr&Good'` (default), `'Ftr&AP&Good'`).
- `sifted.rho`/`pval` are the correlation-filter thresholds (not the
  tutorial paper's 0.3); `K` is the cluster count for the GA-scored
  clustering stage (see SKILL.md).
- `pilot.ntries` (10) is materially lower than the papers' quoted default
  of 30. `pilot.dims` (2 or 3) replaces the legacy boolean
  `opts.pilot.ISA3D`, still accepted as an alias. `pilot.method`:
  `'standard'` (BFGS/analytic) or `'pls'` (Partial Least Squares,
  code addition beyond either paper). `pilot.alpha` weights the
  performance-reconstruction term relative to the feature term (standard
  method only).
- `cloister.corrThreshold`/`maxFeatures` -- field was renamed from
  `cthres` (still accepted as a legacy alias); `maxFeatures` (default 20)
  is the enumeration guard (see SKILL.md's CLOISTER section).
- `pythia.classifier` is a registry name (`'knn'` default, `'svm'`,
  `'tree'`, `'nb'`, `'linear'`, `'ensemble'`; resolved via
  `ISAgetClassifierFcn`), not always SVM as in the original papers.
  `pythia.tuning`: `'sobol'` (default), `'bayes'`, or `'none'` (requires
  `pythia.params`). `pythia.kFold` was renamed from `cvfolds`.
- `trace.method`: `'trace3'` (default, 2D+3D) or `'legacy'` (pre-refactor
  DBSCAN, 2D only, uses `trace.contra` for contradiction removal --
  ignored by `'trace3'`). **`trace.nn`/`trace.prior` do not exist in the
  current schema** -- they were removed; do not add them expecting an
  effect (see SKILL.md's TRACE section for why).
- `general.parallel`/`.ncores` -- renamed from a separate top-level
  `parallel.flag`/`.ncores` block in older versions. `ncores` (18) is a
  workstation-specific default; set to your actual core count or leave
  `parallel: false`.
- `outputs.fig` (not shown above, default false) -- for 3D projections
  only, also writes a `.fig` file alongside each footprint PNG for
  interactive rotation in MATLAB; open with `ISArecallView` to snap it
  back to the optimised viewpoint.

**Removed fields** (present in older versions, no longer read; a stale
`options.json` containing them is harmless but has no effect):
`opts.pythia.uselibsvm`, `opts.trace.usesim`, `opts.sifted.NTREES`.
`opts.pythia.useknn` is deprecated but still honoured as a fallback.

## The `InstanceSpace` class -- primary API

```matlab
obj = InstanceSpace(rootdir);              % reads options.json if present, else defaults
obj = obj.build();                          % run the full pipeline
obj = obj.explore(testRootDir);              % evaluate a trained model on new data
results = obj.getResults();                  % training results (== obj.model)
testResults = obj.getResults(1);             % first explore() call's results
obj.save();                                  % write rootdir/model.mat (-v7.3)
obj = InstanceSpace.load(rootdir);           % read it back
obj.plot();                                  % regenerate PNGs from obj.model
```

Because `InstanceSpace` is a **value class** (not a handle class), every
mutating method returns the updated object -- always assign the result
back (`obj = obj.build()`), or the change is silently discarded.

Staged usage -- change options between individual stages, and only the
stages that need to re-run do (re-running an earlier stage automatically
invalidates every stage that transitively depends on it, so a
save/explore/output step can never silently mix a freshly-recomputed
stage with stale downstream results):

```matlab
obj = InstanceSpace(rootdir);
obj = obj.build('stages', {'prelim', 'sifted', 'pilot'});
obj.opts.pilot.alpha = 2.0;                  % adjust after inspecting the projection
obj = obj.build('stages', {'pilot'});         % re-runs PILOT only; sifted output is reused
obj = obj.build('stages', {'cloister', 'pythia', 'trace'});
```

Canonical stage order (enforced regardless of how `'stages'` lists them):
`prelim -> sifted -> pilot -> cloister -> pythia -> trace`. `cloister`
and `pythia` both depend only on `pilot`'s output, not on each other, so
either may run first within a single `build()` call; `trace` depends on
`pythia`'s predictions and always runs last among the analysis stages.

`buildIS.m`/`exploreIS.m` remain as thin wrappers around this class for
existing callers (notably the MATILDA web platform); prefer the class
directly for new code. See `help InstanceSpace` for the full method list.

## Running the pipeline (function-call form)

```matlab
model = buildIS('/path/to/rootdir/');   % trailing slash required
```

`buildIS.m` runs, in order: load CSV -> optional `FILTER` subsetting ->
PRELIM -> optional `FILTER` again on the SIFTED-reduced features ->
SIFTED -> PILOT -> CLOISTER -> PYTHIA (unless `opts.pythia.skip`) ->
TRACE -> writes `model` struct -> `scriptcsv`/`scriptpng`/`scriptweb` per
the `outputs` flags. Every stage's full output (including intermediate
matrices) lives in the returned `model` struct under `model.<stage>.*`
(e.g. `model.pilot.Z`, `model.trace.best{i}`).

To project new, unseen instances into an already-built space, use
`exploreIS.m`/`obj.explore()` rather than rerunning `buildIS.m` on the
combined data -- it reuses the trained PILOT/PYTHIA/TRACE parameters
instead of re-fitting them, which is the methodologically correct way to
evaluate generalisation to new instances. `explore()` re-scores existing
TRACE footprints against the new points rather than rebuilding them (see
SKILL.md's TRACE section for what happens with a genuinely new
algorithm).

## Outputs

With `outputs.csv`/`png` true, `scriptcsv.m`/`scriptpng.m` write, per
algorithm: footprint boundary vertices, good/best/beta-hard scatter
plots, PYTHIA's predicted-good-region plots, and a summary CSV table
(area, density, purity for good and best footprints, per algorithm, plus
the space's own baseline row). Column headers adapt automatically to 2D
(`z_1, z_2`) vs 3D (`z_1, z_2, z_3`) projections, and 3D PNGs include a
coloured compass overlay for orientation. `outputs.web` writes
colour-scaled CSVs consumed by MATILDA's web front end (leave `false`
unless building for that platform).

## Migrating an older model or options file

`ISAmigrateModel.m` brings a pre-v0.9.0 `model.mat` forward to the
current field layout: `opts` struct renames/merges (e.g.
`opts.perf.MaxMin` -> `MaxPerf`), the `model.data.bestPerformace` typo
fix, `pythia.svm`/`.knn` -> `.classifiers`, missing `completedStages`
inference, and -- where the original training data (`model.pilot.Z` and
`model.data.Y`/`Ybin`/`Ybest`/`algolabels`) is still present --
**retraining** any classifiers still in the legacy LIBSVM struct format
(no `predict()` method, so they can't just be relabelled) using the
current registry, and **recomputing** any pre-refactor DBSCAN+polyshape
footprints fresh via TRACE3. If training data is missing, migration still
renames fields but leaves legacy classifiers in place with a warning;
such a model can still be evaluated via `explore()`/PYTHIA eval mode
(which dispatches to `svmpredict` when it detects a legacy struct instead
of a MATLAB classifier object), but this requires the LIBSVM MEX-files,
which are **not bundled with this repository** (removed as precompiled
binaries with no corresponding source in the tree -- see #29). `PYTHIA`
raises a clear `ISA:PYTHIA:noLibsvm` error naming the missing dependency
if you hit this path without them; obtain LIBSVM from
[the official project](https://www.csie.ntu.edu.tw/~cjlin/libsvm/) if you
actually need to evaluate such a model rather than retraining it.

```matlab
ISAmigrateModel(rootdir);   % migrates rootdir/model.mat in place, with automatic backup
```

## Testing/examples in this repo (useful as reference)

- `example.m` -- minimal, single run of `buildIS`+`exploreIS` on the
  bundled `test/data/` reference dataset with sensible defaults; a
  handful of commonly-adjusted settings (classifier, tuning strategy,
  projection dimensionality, feature selection on/off) exposed as plain
  variables near the top. Outputs land in `test/data/example/`. Start
  here.
- `test_integration.m` -- exhaustive option-coverage regression suite:
  every classifier, tuning strategy, 2D/3D, PLS, viewpoint groups, staged
  `build()`/`explore()`/save-load round-trips through the class API, and
  the full `ISAmigrateModel` legacy-migration table, each in its own
  `test/data/<case_name>/` subdirectory. A good reference for how a
  specific option is meant to be used, but not the place to start reading
  -- edit `defaultOpts()` for shared settings or a specific case's
  `override` function for that case only. **`options.json` at the repo
  root is a generated artifact of these scripts, not a source file** --
  hand-editing it has no lasting effect, since the next run overwrites it
  from the `opts` struct built in MATLAB.
- `liveDemoIS.m` -- stage-by-stage interactive walkthrough through the
  `InstanceSpace` class, `%%`-sectioned for MATLAB's Live Editor (open
  and run cell-by-cell, or "Save As -> MATLAB Live Code File" for a
  `.mlx`). Good for demoing individual stage outputs (e.g. plotting
  `obj.model.pilot.Z` right after the PILOT stage) without writing a
  script from scratch.

## Related, not-yet-integrated repository

`CSimpson4224/InstanceSpace3D` is the original standalone implementation
accompanying the ISA3 paper (Simpson et al. 2025); its 3D/TRACE3 logic
appears to have been substantially merged into this main
`andremun/InstanceSpace` repo (the `opts.pilot.dims=3` branches in
`PILOT.m`/`PILOTviewpoint.m` and the `trace3` default method in
`TRACE.m`), but if a specific ISA3 paper figure needs exact reproduction,
prefer the paper's own repository and diff against this one rather than
assuming full parity.

## Licence

**PolyForm Noncommercial License 1.0.0.** Use, modification, and
redistribution of this toolkit are permitted for non-commercial research
and educational purposes only; commercial use requires prior written
permission from the copyright holders. Every source file carries an
`SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0` header;
the full text is at the repository root (`LICENSE`), with citation
metadata in `CITATION.cff`. This changed from GPL v3 as of v0.9.0 -- a
pre-v0.9.0 copy obtained under GPL v3 remains under GPL v3 as licensed,
but the toolkit as it exists in this repository today, and this
reference's description of it, are covered by PolyForm Noncommercial.

## Acknowledgements

Funding for the development of this toolkit was provided by:

- The Australian Research Council, through the Australian Laureate
  Fellowship FL140100012.
- The University of Melbourne, through grant 2025DYA013.
- The Australian Research Council, through the ARC Industrial
  Transformation Training Centre in Optimisation Technologies, Integrated
  Methodologies, and Applications (OPTIMA; grant No. IC200100009).
