# PRELIM

Label good performance and normalise features and performance

## Syntax

```
[X,Y,out] = PRELIM(X,Y,opts)
[X,Y,out] = PRELIM(X,Y,opts,trainedPrelim)
```

## Description

`[X,Y,out] = PRELIM(X,Y,opts)` fits the preprocessing to the training data and applies it. PRELIM labels each algorithm as *good* or *not good* on each instance, finds the best algorithm per instance, bounds feature outliers, and applies Box-Cox and Z-score transforms to the features and the performance. The fitted bounds and transform parameters are returned in `out`.

`[X,Y,out] = PRELIM(X,Y,opts,trainedPrelim)` applies the bounds and transforms of an earlier training call, `trainedPrelim`, to new data instead of fitting new ones. The good/best labelling is the same code in both modes, so training and evaluation label ties the same way.

## Examples

### Preprocess a performance matrix

Misclassification error is a cost (`MaxPerf = false`). Call an algorithm good when its error is below 0.20.

```matlab
T = readtable('test/data/metadata.csv');
X = T{:, startsWith(T.Properties.VariableNames, 'feature_')};
Y = T{:, startsWith(T.Properties.VariableNames, 'algo_')};

opts = struct('MaxPerf', false, 'AbsPerf', true, 'epsilon', 0.20, ...
              'betaThreshold', 0.55, 'auto', true, 'bound', true, ...
              'norm', true, 'iqrMultiplier', 5);
[Xn, Yn, out] = PRELIM(X, Y, opts);

fprintf('%.0f%% of instances are beta-easy\n', 100*mean(out.beta))
histogram(out.numGoodAlgos)
```

### Apply a trained preprocessing to new instances

```matlab
T2 = readtable('test/data/metadata_test.csv');
X2 = T2{:, startsWith(T2.Properties.VariableNames, 'feature_')};
Y2 = T2{:, startsWith(T2.Properties.VariableNames, 'algo_')};
[X2n, Y2n, out2] = PRELIM(X2, Y2, opts, out);
```

## Input Arguments

### `X` — Feature matrix

*numeric matrix*

One row per instance, one column per feature. May contain `NaN`.

### `Y` — Performance matrix

*numeric matrix*

One row per instance, one column per algorithm. May contain `NaN`; a missing value never counts as good.

### `opts` — Preprocessing options

*structure*

When PRELIM is called by `InstanceSpace`, this structure is assembled from `opts.perf`, `opts.auto`, `opts.bound`, `opts.norm` and `opts.prelim`.

#### `opts.MaxPerf` — Performance direction

*`false` (default) | `true`*

`true` if larger performance values are better (for example accuracy). `false` for a cost, such as run time or error. Set from `opts.perf.MaxPerf`.

#### `opts.AbsPerf` — Threshold type

*`false` (default) | `true`*

`true`: an algorithm is good when its performance is better than `epsilon`. `false`: an algorithm is good when its performance is within a fraction `epsilon` of the best algorithm on that instance. Set from `opts.perf.AbsPerf`.

#### `opts.epsilon` — Good-performance threshold

*`0.05` (default) | scalar*

Absolute threshold, or relative tolerance, depending on `AbsPerf`. Set from `opts.perf.epsilon`.

#### `opts.betaThreshold` — Easy-instance fraction

*`0.55` (default) | scalar in [0, 1]*

An instance is *beta-easy* when more than this fraction of the algorithms are good on it. Set from `opts.perf.betaThreshold`.

#### `opts.auto` — Run the preprocessing

*`true` (default) | `false`*

`false` skips the outlier bounds and the transforms; only the labelling runs. Set from `opts.auto.preproc`.

#### `opts.bound` — Bound outliers

*`true` (default) | `false`*

Clip each feature to `median ± iqrMultiplier*IQR`. Set from `opts.bound.flag`.

#### `opts.norm` — Normalise

*`true` (default) | `false`*

Apply Box-Cox and then Z-score transforms to `X` and `Y`. Set from `opts.norm.flag`.

#### `opts.iqrMultiplier` — Outlier bound width

*`5` (default) | positive scalar*

Set from `opts.prelim.iqrMultiplier`.

### `trainedPrelim` — Trained preprocessing

*structure*

The `out` of an earlier training-mode call. Its presence selects evaluation mode.

## Output Arguments

### `X` — Preprocessed features

*numeric matrix*

### `Y` — Preprocessed performance

*numeric matrix*

With `AbsPerf = false`, `Y` is first converted to performance relative to the best algorithm, then normalised.

### `out` — Labels and fitted parameters

*structure*

#### `out.Ybin` — Good-performance labels

`ninst`-by-`nalgos` logical matrix.

#### `out.Ybest` — Best performance

Best raw performance on each instance.

#### `out.P` — Best algorithm

Index of the best algorithm on each instance. Ties are broken at random.

#### `out.numGoodAlgos` — Number of good algorithms

Per instance.

#### `out.beta` — Beta-easy instances

Logical vector; true when more than `betaThreshold` of the algorithms are good.

#### `out.hibound`, `out.lobound`, `out.medval`, `out.iqrange` — Outlier bounds

Training mode only.

#### `out.minX`, `out.lambdaX`, `out.muX`, `out.sigmaX` — Feature transform

Training mode only. Shift, Box-Cox parameter, mean and standard deviation per feature.

#### `out.minY`, `out.lambdaY`, `out.muY`, `out.sigmaY` — Performance transform

Training mode only. `minY` is one shift for the whole matrix; the others are per algorithm.

## Algorithms

1. **Labelling.** Missing performance values are treated as the worst possible value. With `MaxPerf` and `AbsPerf` both false, algorithm *j* is good on instance *i* when `Y(i,j)/Ybest(i) - 1 <= epsilon`.
2. **Bounding.** Each feature is clipped to `median ± iqrMultiplier*IQR`.
3. **Normalisation.** Each feature is shifted to a minimum of 1, and the whole performance matrix to a small positive minimum. Each column is then Box-Cox transformed with its own fitted λ and standardised to zero mean and unit variance.

In evaluation mode, PRELIM warns (`ISA:InstanceSpace:outOfDistribution`) when more than 5% of the new instances fall outside the training bounds.

## Version History

### v0.9.1 — Evaluation mode

The fourth argument, `trainedPrelim`, applies a trained preprocessing to new data. `InstanceSpace.explore` uses it, so training and evaluation share one implementation.

### v0.9.0 — Outlier bound option

`opts.prelim.iqrMultiplier` replaces the fixed outlier-bound width.

## References

- Smith-Miles, K. & Muñoz, M.A. (2023). Instance Space Analysis for Algorithm Testing. *ACM Computing Surveys*, 55(12), Article 255. <https://doi.org/10.1145/3572895>

## See Also

`INIT` | `FILTER` | `SIFTED` | `InstanceSpace`
