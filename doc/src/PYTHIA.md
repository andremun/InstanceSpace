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

#### `out.Pr0hat` — Predicted bad-class scores

Bad-class score for each algorithm on each instance. Consult the corresponding `out.scoreType` entry before interpreting a column: only `probability` denotes the probability that the algorithm is *not* good. `decision-score` and `class-score` are classifier scores and need not lie in [0,1]; `unknown` has no declared score semantics, and `unavailable` marks a placeholder rather than a usable score.

#### `out.Ysub`, `out.Pr0sub` — Cross-validated predictions

Training mode only. `Ysub` contains predicted labels; `Pr0sub` contains bad-class scores. `out.scoreTypeCV` describes each score column. A `mixed` column combines score types across folds; `out.Pr0subIsProbability` identifies which individual entries can be interpreted as probabilities.

#### `out.scoreType`, `out.scoreTypeCV` — Score semantics

Cell arrays with one entry per algorithm. Values are `probability`, `decision-score`, `class-score`, `unknown`, or `unavailable`. `scoreTypeCV` can also be `mixed`. Use this metadata when consuming `Pr0hat` or `Pr0sub`, including models loaded from older toolkit versions.

#### `out.accuracy`, `out.precision`, `out.recall`, `out.cvcmat` — Classifier performance

Per algorithm. Training mode: from cross-validation. Evaluation mode: on the new instances with observed performance; `NaN` for an algorithm the new data does not cover or that has no trained classifier.

#### `out.param1`, `out.param2` — Selected hyperparameters

Per algorithm. Training mode only.

#### `out.selection0`, `out.selection1` — Selected algorithm

Index of the selected algorithm per instance. `selection0` is 0 where no algorithm is predicted good; `selection1` uses the algorithm with the most good instances there instead.

#### `out.mu`, `out.sigma` — Standardisation of Z

#### `out.summary` — Summary table

Cell array with one row per algorithm plus rows for the *Oracle* (always the best algorithm) and the *Selector*. Columns: mean and standard deviation of performance on all instances, probability of good performance (over the instances with observed performance), mean and standard deviation on the instances where the algorithm is selected, CV accuracy, precision and recall during training, test metrics during exploration, and (training mode) the hyperparameters. Cells with no data, such as the accuracy of an algorithm without a classifier, are empty (`[]`). Written to `classifier_table.csv` by `scriptcsv`.

## GNU Octave

The validated Octave classifier is KNN with none or Sobol tuning, weighted/unweighted fitting, CV, probabilities and held-out evaluation. Other classifiers and Bayesian tuning fail explicitly. The tested environment is Octave 11.3, Statistics 2.0.0 and Datatypes 1.5.0. Actual Sobol candidate matrices are stored in tuningCandidates; equal seeds do not imply the same MATLAB scramble. Versioned model archives preserve native KNN models through the package serialization API.

## Version History

### v0.9.2 — Evaluation corrections and Octave KNN

Added Octave KNN with none/Sobol tuning, stored candidate matrices and classifier archive support.

Evaluation uses the saved training fallback algorithm and precision weights. Legacy models without these fields use the first trained algorithm as fallback and equal voting weights. Test outcomes never determine recommendations.

Training summaries use out-of-fold predictions for both algorithm and selector rows. `selection0CV` and `selection1CV` retain these selections. `Yhat`, `selection0`, and `selection1` remain fitted-data outputs for footprints and training plots. These CV metrics condition on the fitted preprocessing, projection, and selected hyperparameters. They are not an unbiased cross-validation estimate of the complete ISA pipeline. Exploration summaries label their metrics as `Test_model_*`.

Selector recall is the fraction of instances with an observed good algorithm on which the non-fallback selection is good. Successful selections are not also counted as missed opportunities when other algorithms are good.

Cost-sensitive weights are `abs(Y-Ybest)` per instance. Zero regrets use the smallest positive regret in the training matrix to keep every observed example trainable. If all regrets are zero, weights are uniform.

`scoreType` describes each `Pr0hat` column and `scoreTypeCV` describes `Pr0sub`. Values are `probability`, `decision-score`, `class-score`, or `unavailable` (`unknown` for old classifiers without metadata). Scores are mapped through classifier class names. SVM folds and the final model retain posterior calibration. A failed calibration is labelled as decision scores. `Pr0subIsProbability` identifies calibrated or probabilistic CV scores per prediction, and `scoreTypeCV` is `mixed` when folds use different score types. Failed CV candidates cannot produce a successful model: all failed candidates or an invalid selected CV result raise an error. A single-class training fold predicts its observed class.

Evaluation mode scores only the instances with observed performance for each algorithm, and reports `NaN` for a trained algorithm that the new data does not cover. The Oracle's probability of good performance is the fraction of instances on which any algorithm is good, instead of a fixed 1. With `opts.classifier = 'svm'`, the posterior probabilities (`Pr0hat`/`Pr0sub`) are now reproducible for a fixed `opts.seed` — the classifier is reseeded immediately before fitting its sigmoid calibration, rather than inheriting whatever random state `fitcsvm`'s solver left behind. This changes `Pr0hat`/`Pr0sub` for existing SVM runs, but not the thresholded predictions (`Yhat`/`Ysub`).

### v0.9.0 — Classifier registry

`opts.classifier` selects any registered classifier, with Sobol or Bayesian tuning. Replaces the LIBSVM-based SVM and `PYTHIA2`.

## References

- Smith-Miles, K. & Muñoz, M.A. (2023). Instance Space Analysis for Algorithm Testing. *ACM Computing Surveys*, 55(12), Article 255. <https://doi.org/10.1145/3572895>

## See Also

`ISAgetClassifierFcn` | `TRACE` | `InstanceSpace` | [Deprecated Functions](Deprecated.html)
