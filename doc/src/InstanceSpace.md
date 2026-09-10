# InstanceSpace

## Overview

`InstanceSpace` is a MATLAB value class that wraps the ISA pipeline (PRELIM, SIFTED, PILOT, CLOISTER, PYTHIA, TRACE). This class provides a structured way to manage and analyze data through a series of predefined stages. Because `InstanceSpace` is a value class (not a handle class), every method that changes state returns the updated object; the caller must reassign, e.g., `obj = obj.build()`.

`buildIS.m` and `exploreIS.m` are backward-compatibility wrappers around this class; new code is encouraged to use the class directly.

## Construction and Persistence

### Constructor

The constructor for `InstanceSpace` is designed to accept two main arguments: `rootdir` and `opts`. The constructor reads `metadata.csv` and `options.json` to initialize the object, validates the options using `ISAvalidateOpts`, and fills in defaults using `ISAdefaults`.

- **Syntax:** `obj = InstanceSpace(rootdir)` or `obj = InstanceSpace(rootdir, opts)`.
- **Persistence:** The constructor does not run any computation; it merely initializes the object.

### Saving/Loading

- **`save(obj)`:** Writes `model.mat` (-v7.3) by flattening `obj.model` fields.
- **`load(rootdir)` (Static):** Loads `model.mat`, migrates legacy names via `ISAmigrateModel`, and fills missing defaults.

## Properties

The following properties are publicly accessible:

- `rootdir` (1,:) char -- the directory this object was constructed from
- `opts` (1,1) struct -- the full pipeline options struct
- `model` (1,1) struct -- the trained model, populated stage by stage as `build()` runs
- `testDirs` (1,:) cell -- one entry per `explore()` call, the test root directory used
- `testResults` (1,:) cell -- one entry per `explore()` call, that call's result struct
- `completedStages` (1,:) cell -- which of `prelim`, `sifted`, `pilot`, `cloister`, `pythia`, `trace` have run and are still valid; re-running an earlier stage invalidates subsequent stages

## Core Methods

### `build(obj, varargin)`

The `build()` method runs the full pipeline or specific stages.

- **Syntax:** `obj = obj.build()` runs every stage. `obj.build('stages', {'pilot',...})` runs only named stages.
- **Callback:** `obj.build('onStage', @(stageName, model) ...)` invokes a callback once per completed stage with that stage's name and `obj.model` at that point.
- **Persistence:** Outputs are persisted and written (CSV/web/PNG) only once every stage has completed.

### `explore(obj, testRootDir, varargin)`

The `explore()` method evaluates the trained model on new data using frozen training options. It excludes the `cloister` stage.

- **Syntax:** `obj = obj.explore(testRootDir, 'onStage', @(stageName, out) ...)` fires once per conceptual stage.
- **Exclusion:** The `cloister` stage is never recomputed at explore time.

### `getResults(obj, idx)`

- `results = obj.getResults()`: Returns training results (`obj.model`).
- `results = obj.getResults(idx)`: Returns the `idx`-th explore result (`obj.testResults{idx}`).

### `plot(obj, varargin)`

Interactive convenience wrapper around `scriptfcn.m`'s drawing helpers, plotting to the current figure. Requires a trained model.

- `plot(obj, viewName)` -- basic form.
- `plot(obj, viewName, algoIdx)` -- with an algorithm index, required for the `'good'`/`'footprint'` views.

**View options:**

| `viewName` | Notes |
|---|---|
| `'sources'` | Requires `model.data.S` |
| `'portfolio'` | Uses `drawPortfolioSelections` |
| `'good'` | Requires `algoIdx` (1-based into `model.data.algolabels`) |
| `'footprint'` | Requires `algoIdx` (1-based into `model.data.algolabels`) |
| `'boundary'` | Requires `model.cloist`; 2D only |

## Usage

### Basic usage

```matlab
obj = InstanceSpace(rootdir);
obj = obj.build();
obj = obj.explore(testRootDir);
results = obj.getResults();
testResults = obj.getResults(1);
```

### Staged usage, with option changes between stages

```matlab
obj = InstanceSpace(rootdir);
obj.opts.pilot.dims = 3;
obj = obj.build('stages', {'prelim','sifted','pilot'});
figure; scatter(obj.model.pilot.Z(:,1), obj.model.pilot.Z(:,2));
obj.opts.pilot.alpha = 2.0;
obj = obj.build('stages', {'pilot'});
obj.opts.pythia.tuning = 'bayes';
obj = obj.build('stages', {'cloister','pythia','trace'});
```

Re-running an earlier stage invalidates every later stage's output.

### Callback-based inspection

```matlab
inspect = @(stageName, model) fprintf('%s done\n', stageName);
obj = obj.build('onStage', inspect);
obj = obj.explore(testRootDir, 'onStage', inspect);
```

## References

- Smith-Miles, K. & Munoz, M.A. (2023). Instance Space Analysis for Algorithm Testing. ACM Computing Surveys, 55(12), Article 255. <https://doi.org/10.1145/3572895>

## See Also

[buildIS](buildIS.html) | [exploreIS](exploreIS.html)
