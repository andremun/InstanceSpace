# PYTHIA

Train classifiers that predict which algorithms perform well

<!-- opts: pythia -->

## Syntax

```
out = PYTHIA(Z,Y,Ybin,Ybest,algolabels,opts)
out = PYTHIA(Z,Y,Ybin,Ybest,algolabels,opts,trainedModel)
```

## Description

`out = PYTHIA(Z,Y,Ybin,Ybest,algolabels,opts)` trains one binary classifier per algorithm that predicts, from an instance's coordinates `Z`, whether the algorithm is good on it (Performance Yielding Heuristic Identification of Algorithms). Hyperparameters are tuned by *k*-fold cross-validation. PYTHIA then combines the predictions into an algorithm selector: for each instance it picks the algorithm predicted to be good that has the highest cross-validated precision.

`out = PYTHIA(Z,Y,Ybin,Ybest,algolabels,opts,trainedModel)` applies the classifiers of an earlier training call to new instances and reports how accurate they are there. Classifiers are not retrained.

The classifier type is set by `opts.classifier`; see `ISAgetClassifierFcn` for the registry and the tuned hyperparameters.

## Examples

### Train k-NN classifiers on a 2D instance space

```matlab
rootdir = 'test/data/example/';
if ~isfolder(rootdir), mkdir(rootdir); end
copyfile('test/data/metadata.csv', rootdir);
opts.perf = struct('MaxPerf', false, 'AbsPerf', true, 'epsilon', 0.20);
obj = InstanceSpace(rootdir, opts);
obj = obj.build('stages', {'prelim', 'sifted', 'pilot'});
d = obj.model.data;

out = PYTHIA(obj.model.pilot.Z, d.Yraw, d.Ybin, d.Ybest, d.algolabels, obj.opts.pythia);
disp(out.summary)
```

### Use an SVM with Bayesian optimisation

```matlab
pythiaOpts = obj.opts.pythia;
pythiaOpts.classifier = 'svm';
pythiaOpts.tuning = 'bayes';
pythiaOpts.nTuningIter = 30;
out = PYTHIA(obj.model.pilot.Z, d.Yraw, d.Ybin, d.Ybest, d.algolabels, pythiaOpts);
```

### Evaluate trained classifiers on new instances

`InstanceSpace.explore` does this for you. The explicit call is:

```matlab
copyfile('test/data/metadata_test.csv', rootdir);
obj = obj.build();                  % train every stage
obj = obj.explore(rootdir);         % reads rootdir/metadata_test.csv
res = obj.getResults(1);
res.pythia.summary
```

## Input Arguments

### `Z` — Instance coordinates

*numeric matrix*

`model.pilot.Z`. PYTHIA standardises it internally and stores the mean and standard deviation.

### `Y` — Raw performance matrix

*numeric matrix*

`model.data.Yraw`. Used for the performance columns of the summary and, with `opts.useweights`, for the class weights. `NaN` marks a missing value.

### `Ybin` — Good-performance labels

*logical matrix*

The classification targets.

### `Ybest` — Best performance per instance

*numeric vector*

### `algolabels` — Algorithm names

*cell array of character vectors*

### `opts` — Classifier options

*structure*

Normally `obj.opts.pythia`.

#### `opts.classifier` — Classifier type

*`'knn'` (default) | `'svm'` | `'tree'` | `'nb'` | `'linear'` | `'ensemble'`*

#### `opts.tuning` — Hyperparameter search

*`'sobol'` (default) | `'bayes'` | `'none'`*

`'sobol'`: a scrambled Sobol sequence of `nTuningIter` candidates. `'bayes'`: `bayesopt` with the same budget. `'none'`: use `opts.params`.

#### `opts.nTuningIter` — Search budget

*`20` (default) | positive integer*

#### `opts.kFold` — Cross-validation folds

*`5` (default) | positive integer*

#### `opts.params` — Fixed hyperparameters

*`[]` (default) | `nalgos`-by-`nparams` matrix*

Required when `opts.tuning` is `'none'`; one row per algorithm.

#### `opts.useweights` — Cost-sensitive training

*`false` (default) | `true`*

Weight each instance by how far its performance is from the mean performance.

#### `opts.ispolykrnl` — Polynomial kernel

*`false` (default) | `true`*

SVM only. `false` uses a Gaussian kernel.

#### `opts.ensembleMethod` — Ensemble method

*`'Bag'` (default) | any `fitcensemble` method*

Only for `'ensemble'`.

#### `opts.skip` — Skip training

*`false` (default) | `true`*

Return empty predictions without training. `TRACE` then uses `Ybin` directly.

#### `opts.seed` — Random seed

*`opts.general.seed` (default) | integer*

#### `opts.verbose` — Report progress

*`opts.general.verbose` (default) | logical*

### `trainedModel` — Trained classifiers

*structure*

The `out` of an earlier training call (`model.pythia`). Its presence selects evaluation mode.

## Output Arguments

### `out` — Classifiers, predictions and summary

*structure*

#### `out.classifiers` — Trained classifiers

Cell array with one `ClassificationModel` per algorithm.

#### `out.Yhat` — Predicted labels

`ninst`-by-`nalgos` logical matrix. In training mode these are the predictions of the final model on the training data.

#### `out.Pr0hat` — Predicted probabilities

Probability that each algorithm is *not* good on each instance.

#### `out.Ysub`, `out.Pr0sub` — Cross-validated predictions

Training mode only.

#### `out.accuracy`, `out.precision`, `out.recall`, `out.cvcmat` — Classifier performance

Per algorithm. Training mode: from cross-validation. Evaluation mode: on the new instances with observed performance; `NaN` for an algorithm the new data does not cover or that has no trained classifier.

#### `out.param1`, `out.param2` — Selected hyperparameters

Per algorithm. Training mode only.

#### `out.selection0`, `out.selection1` — Selected algorithm

Index of the selected algorithm per instance. `selection0` is 0 where no algorithm is predicted good; `selection1` uses the algorithm with the most good instances there instead.

#### `out.mu`, `out.sigma` — Standardisation of Z

#### `out.summary` — Summary table

Cell array with one row per algorithm plus rows for the *Oracle* (always the best algorithm) and the *Selector*. Columns: mean and standard deviation of performance on all instances, probability of good performance (over the instances with observed performance), mean and standard deviation on the instances where the algorithm is selected, cross-validation accuracy, precision and recall, and (training mode) the hyperparameters. Cells with no data, such as the accuracy of an algorithm without a classifier, are empty (`[]`). Written to `classifier_table.csv` by `scriptcsv`.

## Version History

### v0.9.2 — Evaluation, Oracle, and SVM reproducibility fixes

Evaluation mode scores only the instances with observed performance for each algorithm, and reports `NaN` for a trained algorithm that the new data does not cover. The Oracle's probability of good performance is the fraction of instances on which any algorithm is good, instead of a fixed 1. With `opts.classifier = 'svm'`, the posterior probabilities (`Pr0hat`/`Pr0sub`) are now reproducible for a fixed `opts.seed` — the classifier is reseeded immediately before fitting its sigmoid calibration, rather than inheriting whatever random state `fitcsvm`'s solver left behind. This changes `Pr0hat`/`Pr0sub` for existing SVM runs, but not the thresholded predictions (`Yhat`/`Ysub`).

### v0.9.0 — Classifier registry

`opts.classifier` selects any registered classifier, with Sobol or Bayesian tuning. Replaces the LIBSVM-based SVM and `PYTHIA2`.

## References

- Smith-Miles, K. & Muñoz, M.A. (2023). Instance Space Analysis for Algorithm Testing. *ACM Computing Surveys*, 55(12), Article 255. <https://doi.org/10.1145/3572895>

## See Also

`ISAgetClassifierFcn` | `TRACE` | `InstanceSpace` | [Deprecated Functions](Deprecated.html)
