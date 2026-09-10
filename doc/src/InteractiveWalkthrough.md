# InteractiveWalkthrough

## Introduction

An interactive, stage-by-stage walkthrough of the ISA pipeline using the `InstanceSpace` class. Adapted from `liveDemoIS.m`, a `%%`-sectioned file designed to be opened in MATLAB's Live Editor (distinct from `example.m`, the Getting Started page's source).

## Installation Requirements

To run this walkthrough, ensure you have the following:
- MATLAB R2025a or later.
- The following toolboxes installed:
  - Global Optimization
  - Parallel Computing
  - Optimization
  - Statistics and Machine Learning
  - Financial (for `boxcox()`).
- Communications Toolbox is not required.
- LIBSVM support is deprecated; new runs use the native classifier registry.

Before starting, run `startup.m` once per session to ensure all necessary paths are added. Alternatively, you can construct the `InstanceSpace` class to add paths automatically.

## Setup

The walkthrough uses a reference dataset from Munoz et al. (2018), which consists of 212 training instances and 23 test instances, with 10 features and 10 algorithms. The dataset's performance scores are measured by misclassification error, where a lower score indicates better performance.

### Setup Process

1. Copy `metadata.csv` and `metadata_test.csv` into a root directory.
2. Set specific `opts.perf` values:
   - `MaxPerf=false`
   - `AbsPerf=true`
   - `epsilon=0.20`
3. Construct the `InstanceSpace` object: `obj = InstanceSpace(rootdir, opts)`.

After construction, `obj.completedStages` is empty, indicating no stages have been completed yet.

## PRELIM Stage

### Process Description

The PRELIM stage involves:
- Removing instances/features with too many missing values.
- Optionally bounding outliers (median +/- `opts.prelim.iqrMultiplier` x IQR).
- Normalizing every feature/performance column using Box-Cox and Z-score.
- Determining which algorithm is "good" on which instance based on `opts.perf.*`.

### Execution

Run `obj.build('stages', {'prelim'})` to execute the PRELIM stage. The demo output includes instance/feature/algorithm counts and the fraction of "easy" instances (`mean(obj.model.data.beta)`).

## SIFTED Stage

### Feature Selection Process

Feature selection involves:
- A correlation filter to identify the most predictive features.
- A genetic algorithm (if more than a handful of features remain) to pick representative features from each correlation cluster.

### Execution

If `opts.sifted.flag=false`, skip this stage, using every feature. Run `obj.build('stages', {'sifted'})` to execute the SIFTED stage. The demo output compares feature lists before and after selection (`obj.model.featsel.labels` vs `obj.model.data.featlabels`).

## PILOT Stage

### Linear Projection Explanation

The PILOT stage involves finding a linear projection `Z = X*A'` that reconstructs both features and performance as closely as possible. The projection is used to explore the instance space in 2D or 3D dimensions.

### Options

- `opts.pilot.dims=3` for a 3D projection.
- `opts.pilot.method='pls'` for using Partial Least Squares instead of the default BFGS/analytic method.

### Execution

Run `obj.build('stages', {'pilot'})` to execute the PILOT stage. The demo output includes a scatter plot of `obj.model.pilot.Z(:,1)` vs `Z(:,2)`, colored by `obj.model.data.numGoodAlgos`.

## CLOISTER

**Purpose:** Estimating the reachable boundary of the instance space based on feature correlations.

**Usage:** Running the stage via `obj = obj.build('stages', {'cloister'})`.

**Visualization:** Using `obj.plot('boundary')` to view the boundary.

## PYTHIA

**Purpose:** Training binary classifiers per algorithm to predict performance over instance space `Z`.

**Configuration:** Details on `opts.pythia.classifier` (default 'knn'; others: 'svm', 'tree', 'nb', 'linear', 'ensemble') and `opts.pythia.tuning` (default scrambled Sobol; others: 'bayes' or pre-supplied params).

**Execution & Output:** Running via `obj = obj.build('stages', {'pythia'})`. Display of `obj.model.pythia.summary` (accuracy/precision/recall) and scatter-plots of `Z` colored by `obj.model.pythia.selection0`.

## TRACE

**Purpose:** Building footprints (regions of expected high performance) for individual algorithms and the portfolio.

**Configuration:** Mention of `opts.trace.method='legacy'` for the DBSCAN + alpha-shape algorithm (2D only).

**Execution & Output:** Running via `obj = obj.build('stages', {'trace'})`. Display of `obj.model.trace.summary` (area/density/purity).

## Post-processing

**Saving:** `obj.save()` writes `model.mat` (HDF5-compatible, -v7.3) containing all stage outputs and geometry objects.

**Side Effects:** Note that CSV and PNG outputs are written during `build()`, not `save()`.

## Exploring the model

**Purpose:** Projecting new instances (from `metadata_test.csv`) into the fitted space for evaluation without retraining.

**Execution:** Running `obj = obj.explore(rootdir)` and retrieving results via `testResults = obj.getResults(1)`.

**Visualization:** Display of `testResults.pythia.summary`, `testResults.trace.summary`, and scatter-plots of `testResults.pilot.Z` colored by `testResults.pythia.selection0`.

## Additional resources

**Documentation/Testing:** `README.md` (options/layout) and `test_integration.m` (regression suite/option reference).

**Web Platform:** MATILDA (https://matilda.unimelb.edu.au).

## See Also

[GettingStarted](GettingStarted.html) | [Landing](Landing.html)
