# PYTHIA
Unified binary classifier for instance-space algorithm selection.
## Syntax
```
# Training mode
out = PYTHIA(Z, Y, Ybin, Ybest, algolabels, opts)

# Evaluation mode
out = PYTHIA(Z, Y, Ybin, Ybest, algolabels, opts, trainedModel)
```
## Description
Trains one binary good/not-good classifier per algorithm (registry resolved via ISAgetClassifierFcn), tuning hyperparameters per opts.tuning. Applies trainedModel's already-fitted classifiers to new data instead of training fresh ones -- trainedModel is a prior training-mode call's out struct.
## Input Arguments
| Argument | Description |
|---|---|
| Z | (ninst x ndim) projected instance coordinates |
| Y | (ninst x nalgos) raw performance matrix |
| Ybin | (ninst x nalgos) logical good-performance labels |
| Ybest | (ninst x 1) best algorithm's raw performance per instance |
| algolabels | (1 x nalgos) cell array of algorithm name strings |
| opts | struct with fields: |
| &nbsp;&nbsp;`opts.classifier` | char registry name (default 'knn'); see ISAgetClassifierFcn |
| &nbsp;&nbsp;`opts.tuning` | char 'sobol' (default) -- scrambled Sobol quasi-random search, opts.nTuningIter evals; 'bayes' -- MATLAB bayesopt (Gaussian process surrogate), same evals/k-fold CV as 'sobol'; 'none' -- use pre-supplied opts.params directly, skip tuning |
| &nbsp;&nbsp;`opts.nTuningIter` | int Sobol/Bayes evaluation budget (default 20) |
| &nbsp;&nbsp;`opts.kFold` | int cross-validation folds (default 5) |
| &nbsp;&nbsp;`opts.ispolykrnl` | logical SVM only: polynomial vs Gaussian kernel (default false) |
| &nbsp;&nbsp;`opts.useweights` | logical cost-sensitive classification (default false) |
| &nbsp;&nbsp;`opts.params` | double pre-calculated hyperparameters, required when tuning='none' (default []) |
| &nbsp;&nbsp;`opts.skip` | logical bypass training entirely (default false) |
| &nbsp;&nbsp;`opts.ensembleMethod` | char fitcensemble method (default 'Bag') |
| &nbsp;&nbsp;`opts.seed` | double RNG seed, defaults to opts.general.seed |
| &nbsp;&nbsp;`opts.verbose` | logical per-classifier progress output (default true) |
| trainedModel | (optional) a prior training-mode call's out struct; presence selects evaluation mode over training mode |
## Output Arguments
| Argument | Description |
|---|---|
| out | struct with fields: |
| &nbsp;&nbsp;`out.classifiers{nalgos}` | trained classifier objects (or legacy LIBSVM structs) |
| &nbsp;&nbsp;`out.Yhat`, `out.Pr0hat` | predicted good-performance labels and probabilities |
| &nbsp;&nbsp;`out.Ysub`, `out.Pr0sub` | cross-validated predictions (training mode) |
| &nbsp;&nbsp;`out.cvcmat`, `out.accuracy`, `out.precision`, `out.recall` | per-algorithm cross-validation performance |
| &nbsp;&nbsp;`out.param1`, `out.param2`, `out.param2Label` | selected hyperparameters per algorithm |
| &nbsp;&nbsp;`out.selection0`, `out.selection1` | predicted best-algorithm index per instance (precision-weighted, with and without a default fallback) |
| &nbsp;&nbsp;`out.mu`, `out.sigma` | Z-score normalisation parameters applied to Z |
| &nbsp;&nbsp;`out.summary` | cell array display of the above, one row per algorithm |
## References

- Smith-Miles, K. & Munoz, M.A. (2023). Instance Space Analysis for Algorithm Testing. ACM Computing Surveys, 55(12), Article 255. <https://doi.org/10.1145/3572895>

## See Also

[TRACE](TRACE.html) | [InstanceSpace](InstanceSpace.html)
