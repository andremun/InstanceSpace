# PILOT

Project features onto a 2D or 3D instance space (Munoz et al., Mach Learn 2018)

## Syntax

```
out = PILOT(X, Y, featlabels, opts)
```

## Description

finding A/B/C such that Z=X*A' and [X Y] is reconstructed from Z as
closely as possible.

## Input Arguments

| Argument | Description |
|---|---|
| `X` | (`ninst` x `nfeats`) feature matrix |
| `Y` | (`ninst` x `nalgos`) performance matrix |
| `featlabels` | (1 x `nfeats`) cell array of feature name strings |
| `opts` | struct with fields (see [ISAdefaults](ISAdefaults.html) for defaults) |
| &nbsp;&nbsp;`opts.dims` | `int` -- projection dimensionality, 2 or 3 |
| &nbsp;&nbsp;`opts.method` | `char` -- `'standard'` (BFGS/analytic, see below) or `'pls'` (Partial Least Squares via `plsregress`; `opts.alpha` does not apply) |
| &nbsp;&nbsp;`opts.analytic` | `logical` -- use the closed-form eigenvector solution (`method='standard'` only); falls back to numerical if `X` is rank-deficient |
| &nbsp;&nbsp;`opts.ntries` | `int` -- BFGS multi-start restarts (numerical branch only) |
| &nbsp;&nbsp;`opts.alpha` | `double` -- performance-reconstruction cost weight (`method='standard'` only): `min \|\|Xtilde-BrZ\|\|^2 + alpha*\|\|Y-CrZ\|\|^2` |
| &nbsp;&nbsp;`opts.topoWeight` | `double` -- reserved for future use; not wired into the cost function |
| &nbsp;&nbsp;`opts.verbose` | `logical` -- per-trial progress output |
| &nbsp;&nbsp;`opts.precalcAlpha` | (optional) pre-computed full BFGS solution vector, skips optimisation |
| &nbsp;&nbsp;`opts.X0` | (optional) user-supplied BFGS starting points |

## Output Arguments

| Field | Description |
|---|---|
| `out` | struct with fields: |
| &nbsp;&nbsp;`out.A` | (`dims` x `nfeats`) projection matrix, `Z = X*A'` |
| &nbsp;&nbsp;`out.B`, `out.C` | reconstruction matrices for the feature/performance blocks of `[X Y]` from `Z` |
| &nbsp;&nbsp;`out.Z` | (`ninst` x `dims`) projected instance coordinates |
| &nbsp;&nbsp;`out.error` | sum of squared reconstruction error |
| &nbsp;&nbsp;`out.R2` | per-column R^2 of the reconstruction |
| &nbsp;&nbsp;`out.summary` | cell array display of `A` with feature labels |
| &nbsp;&nbsp;`out.alpha`, `X0`, `eoptim`, `perf` | (numerical branch only) raw BFGS trial results, used to pick the best of `opts.ntries` |

## References

- Smith-Miles, K. & Munoz, M.A. (2023). Instance Space Analysis for Algorithm Testing. *ACM Computing Surveys*, 55(12), Article 255. <https://doi.org/10.1145/3572895>
- Simpson, C., Munoz, M.A., Kandanaarachchi, S. & Campello, R.J.G.B. (2025). ISA3: A 3-dimensional expansion of Instance Space Analysis. *Machine Learning*, 114, 240. <https://doi.org/10.1007/s10994-025-06871-5>
- Munoz, M.A., Villanova, L., Baatar, D. & Smith-Miles, K. (2018). Instance spaces for machine learning classification. *Machine Learning*, 107(1), 109-147. <https://doi.org/10.1007/s10994-017-5629-5>

## See Also

[SIFTED](SIFTED.html) | [PILOTviewpoint](PILOTviewpoint.html) | [InstanceSpace](InstanceSpace.html)
