# PRELIM

Pre-process instance-space data before projection.

## Syntax

### Training mode (3 args)

```
[X, Y, out] = PRELIM(X, Y, opts)
```

### Evaluation mode (4 args)

```
[X, Y, out] = PRELIM(X, Y, opts, trainedPrelim)
```

## Description

Fits outlier bounds and Box-Cox+Z-score normalisation parameters fresh from X/Y and returns them in out.

Applies a previously-fit trainedPrelim's bounds/normalisation parameters to new X/Y instead of re-fitting them -- trainedPrelim is a prior training-mode call's out struct.

## Input Arguments

| Argument | Description |
|---|---|
| `X` | (ninst x nfeats) feature matrix; may contain NaN |
| `Y` | (ninst x nalgos) performance matrix; may contain NaN |
| `opts` | struct with fields: |
| `opts.MaxPerf` | logical true = maximise performance (default false) |
| `opts.AbsPerf` | logical true = absolute threshold (default false) |
| `opts.epsilon` | double good-performance threshold (default 0.05) |
| `opts.betaThreshold` | double easy-instance fraction (default 0.55) |
| `opts.auto` | logical run auto pre-processing (default true) |
| `opts.bound` | logical bound outliers (default true) |
| `opts.norm` | logical Box-Cox + Z normalisation (default true) |
| `opts.iqrMultiplier` | double outlier bound = median +/- N*IQR (default 5) |
| `trainedPrelim` | (optional) a prior training-mode call's out struct; presence selects evaluation mode over training mode |

## Output Arguments

| Argument | Description |
|---|---|
| `X` | pre-processed feature matrix |
| `Y` | pre-processed performance matrix |
| `out` | struct with fields: |
| `out.Ybest`, `out.Ybin`, `out.P`, `out.numGoodAlgos`, `out.beta` | (performance summary) |
| | Training mode only, also: |
| `out.medval`, `out.iqrange`, `out.hibound`, `out.lobound` | (outlier bounds) |
| `out.minX`, `out.lambdaX`, `out.muX`, `out.sigmaX` | (feature normalisation) |
| `out.minY`, `out.lambdaY`, `out.muY`, `out.sigmaY` | (performance normalisation) |

## References

- Smith-Miles, K. & Munoz, M.A. (2023). Instance Space Analysis for Algorithm Testing. ACM Computing Surveys, 55(12), Article 255. <https://doi.org/10.1145/3572895>

## See Also

[INIT](INIT.html) | [FILTER](FILTER.html) | [SIFTED](SIFTED.html) | [InstanceSpace](InstanceSpace.html)
