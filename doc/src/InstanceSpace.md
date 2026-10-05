# InstanceSpace

Build, explore and save an instance space

## Syntax

```
obj = InstanceSpace(rootdir)
obj = InstanceSpace(rootdir,opts)
```

## Description

An `InstanceSpace` object holds the options, the trained model, and the evaluation results of one instance space analysis. It runs the pipeline stages `PRELIM`, `SIFTED`, `PILOT`, `CLOISTER`, `PYTHIA` and `TRACE` in order, lets you change options and re-run individual stages, evaluates the trained model on new instances, and saves and loads the model.

`obj = InstanceSpace(rootdir)` creates an object for the metadata in `rootdir/metadata.csv`. Options are read from `rootdir/options.json` if the file exists; missing options take their defaults. No computation runs until you call `build`.

`obj = InstanceSpace(rootdir,opts)` uses the options structure `opts` instead of `options.json`.

During `build`, algorithms with no good training instances are removed. Preprocessing is then fitted to the retained portfolio. Best-algorithm indices, best performance, beta-easy flags, and saved transforms all refer to that portfolio. With an absolute threshold, removing an algorithm can change the best performance on an instance where no algorithm is good.

`InstanceSpace` is a value class: a method that changes the object returns the changed copy, so assign the result, as in `obj = obj.build()`.

## Examples

### Build and explore an instance space

```matlab
rootdir = 'test/data/example/';
if ~isfolder(rootdir), mkdir(rootdir); end
copyfile('test/data/metadata.csv', rootdir);
copyfile('test/data/metadata_test.csv', rootdir);

opts.perf = struct('MaxPerf', false, 'AbsPerf', true, 'epsilon', 0.20);
obj = InstanceSpace(rootdir, opts);
obj = obj.build();              % every stage; writes model.mat, CSV and PNG files
obj = obj.explore(rootdir);     % evaluates rootdir/metadata_test.csv

model = obj.getResults();       % training results
test  = obj.getResults(1);      % first explore() result
```

### Change an option and re-run one stage

Re-running a stage discards the results of every later stage, which you then re-run.

```matlab
obj = InstanceSpace(rootdir, opts);
obj = obj.build('stages', {'prelim', 'sifted', 'pilot'});
obj.plot('portfolio')

obj.opts.pilot.analytic = true;
obj = obj.build('stages', {'pilot'});
obj = obj.build('stages', {'cloister', 'pythia', 'trace'});
```

### Inspect each stage as it finishes

```matlab
report = @(stage, model) fprintf('%s finished\n', stage);
obj = obj.build('onStage', report);
```

### Save and load a model

```matlab
obj.save();                              % rootdir/model.mat
obj2 = InstanceSpace.load(rootdir);
obj2 = obj2.explore(rootdir);
```

## Input Arguments

### `rootdir` — Data folder

*character vector | string scalar*

Folder containing `metadata.csv`. Output files and `model.mat` are written here. See [Metadata File Format](MetadataFormat.html).

### `opts` — Options

*structure*

Any subset of the fields in the [Options Reference](OptionsReference.html). Values are checked by `ISAvalidateOpts`, and missing fields filled by `ISAdefaults`.

## Properties

### `rootdir` — Data folder

*character vector*

Always ends with a file separator.

### `opts` — Options

*structure*

Complete options. Change fields between `build` calls to re-run stages with new settings. `explore` uses the options stored with the model at training time instead.

### `model` — Trained model

*structure*

One field per completed stage: `data`, `prelim`, `featsel`, `sifted`, `pilot`, `cloist` (CLOISTER), `pythia`, `trace`, and `opts` once every stage has run.

### `completedStages` — Completed stages

*cell array of character vectors*

### `testDirs`, `testResults` — Evaluation results

*cell arrays*

One entry per `explore` call: the folder and the result structure, which has the same fields as `model`.

## Object Functions

| Function | Purpose |
|---|---|
| [`build`](#build) | Run pipeline stages |
| [`explore`](#explore) | Evaluate the trained model on new instances |
| [`getResults`](#getresults) | Return the training or an evaluation result |
| [`plot`](#plot) | Plot a view of the instance space |
| [`save`](#save) | Write `model.mat` |
| [`InstanceSpace.load`](#load) | Create an object from `model.mat` |

### build

```
obj = build(obj)
obj = build(obj,'stages',stages)
obj = build(obj,'onStage',callback)
```

Runs every stage, or only the named `stages` from `{'prelim','sifted','pilot','cloister','pythia','trace'}`, always in pipeline order. A stage whose prerequisites have not completed raises `ISA:InstanceSpace:missingPrereq`. `callback(stageName, model)` is called after each stage.

When all stages have completed, `build` saves `model.mat` and writes the outputs selected in `opts.outputs` with `scriptcsv`, `scriptweb` and `scriptpng`.

### explore

```
obj = explore(obj,testRootDir)
obj = explore(obj,testRootDir,'onStage',callback)
```

Evaluates `testRootDir/metadata_test.csv` with the trained model: the trained preprocessing, feature selection, projection, classifiers and footprints are applied, not refitted. Appends the result to `testResults` and writes the outputs to `testRootDir`. Needs a fully built model.

### getResults

```
results = getResults(obj)
results = getResults(obj,idx)
```

Returns `obj.model`, or the result of the `idx`-th `explore` call.

### plot

```
plot(obj,view)
plot(obj,view,algoIdx)
```

Plots into the current figure. `view` is one of:

| View | Shows | Needs |
|---|---|---|
| `'sources'` | instances coloured by source | a `source` column in the metadata |
| `'portfolio'` | best algorithm per instance | |
| `'good'` | good and bad instances of algorithm `algoIdx` | `algoIdx` |
| `'footprint'` | the footprint of algorithm `algoIdx` over its good and bad instances | `algoIdx` |
| `'boundary'` | the CLOISTER boundary over all instances | the `cloister` stage |

### save

```
save(obj)
```

Writes `obj.model` to `rootdir/model.mat`. MATLAB uses MAT-file version 7.3, one variable per model field. Octave uses a versioned MAT v7 archive with `archiveVersion` and `payload` envelope variables; its trained classifier archives are runtime-specific.

### load

```
obj = InstanceSpace.load(rootdir)
```

Reads `rootdir/model.mat`, migrates a legacy model with `ISAmigrateModel`, fills missing options with `ISAdefaults`, and sets `completedStages` from the stages present. `metadata.csv` is not needed.

## GNU Octave

Octave supports serial build/save/load/explore with the validated KNN/TRACE3 backends. Explicit property setters preserve validation in both engines. Octave uses a versioned MAT v7 archive; native KNN archives are runtime-specific, while MATLAB keeps the existing MAT v7.3 format. CSV and PNG are supported; Qt and a display or virtual display are required for 3D PNG. outputs.fig defaults to false in Octave.

## Version History

### v0.9.2 — Stage validation and Octave workflows

Added serial Octave build/save/load/explore with explicit property validation and runtime-specific model archives. Unsupported output options are rejected before work begins, including in staged builds and explore calls using saved model options. Octave archive schema 2 encloses ordinary structs explicitly so user fields cannot be mistaken for internal type markers; existing schema-1 archives remain readable. Native flattened models with custom `archiveVersion` fields remain native unless their only top-level variables are the reserved archive envelope pair, `archiveVersion` and `payload`.

Each completed stage stores the options used to fit it. A partial rebuild rejects changed options for retained stages. Rebuild the affected stage and its dependents first. Partial builds can be saved, loaded, and resumed. Changing `general.seed` requires rebuilding preprocessing. Stage seeds and verbosity are copied from general defaults at construction and remain explicit stage options afterwards.

SIFTED reruns start from the saved pre-selection data, including its original row subset. Models saved without this artifact must rebuild preprocessing before rerunning SIFTED.

Build reuses an existing pool without replacing or deleting it. A pool opened by build is closed on success or error. `general.parallel=false` makes core stage execution serial even if a pool exists.

### v0.9.1 — Stage callbacks and boundary plot

`build` and `explore` accept `'onStage'`. `plot` adds the `'boundary'` view.

### v0.9.0 — Introduced

Replaces the `buildIS` and `exploreIS` scripts as the main interface.

## References

- Smith-Miles, K. & Muñoz, M.A. (2023). Instance Space Analysis for Algorithm Testing. *ACM Computing Surveys*, 55(12), Article 255. <https://doi.org/10.1145/3572895>

## See Also

`buildIS` | `exploreIS` | [Getting Started](GettingStarted.html) | [Options Reference](OptionsReference.html)
