# CLOISTER
Estimate the empirical boundary of the instance space.
## Syntax
```
out = CLOISTER(X, A, opts)
```
## Description
Enumerates every combination of feature lower/upper bounds (a hypercube's corners), discards combinations that contradict the significant pairwise feature correlations (opts.corrThreshold, opts.pval), and projects the surviving corners through the PILOT projection matrix A to trace a boundary polygon in the instance space. If the feature count exceeds opts.maxFeatures, the corner enumeration is skipped and a plain convex hull of the projected instances is used instead.
## Input Arguments
| Argument | Description |
|---|---|
| `X` | `(ninst x nfeats) feature matrix; may contain sparse NaNs` |
| `A` | `(ndim x nfeats) PILOT projection matrix (model.pilot.A)` |
| `opts` | struct with fields (see ISAdefaults for defaults): |
| `opts.pval` | `double significance level for feature correlations (0.05)` |
| `opts.corrThreshold` | `double |correlation| above which a corner combination contradicting the trend is discarded (0.70)` |
| `opts.maxFeatures` | `int feature-count guard before falling back to a plain convex hull (20)` |
## Output Arguments
| Argument | Description |
|---|---|
| `out` | struct with fields: |
| `out.Zedge` | `boundary polygon vertices using every corner` |
| `out.Zecorr` | `boundary polygon vertices using only correlation-consistent corners (same as Zedge when the correlation threshold rejects nothing, or when maxFeatures triggers the convex-hull fallback)` |

## References

- Smith-Miles, K. & Munoz, M.A. (2023). Instance Space Analysis for Algorithm Testing. ACM Computing Surveys, 55(12), Article 255. <https://doi.org/10.1145/3572895>

## See Also

[PILOT](PILOT.html) | [InstanceSpace](InstanceSpace.html)
