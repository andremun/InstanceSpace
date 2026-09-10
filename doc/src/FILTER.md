# FILTER
Density-based instance subsetting for small-scale experiments.
## Syntax
```
[subsetIndex,isDissimilar,isVISA,unif] = FILTER(X,Y,Ybin,opts)
```
## Description
For every pair of instances closer than opts.mindistance in feature space, the second is marked redundant according to opts.type (see below); the caller (buildIS/InstanceSpace) typically keeps the complement of subsetIndex, i.e. drops the redundant instances.
## Input Arguments
| Argument | Description |
|---|---|
| `X, Y, Ybin` | `(ninst x nfeats)/(ninst x nalgos)/(ninst x nalgos) feature, performance, and good-performance matrices` |
| `opts` | struct with fields: |
| `opts.mindistance` | `double feature-space distance threshold below which two instances are considered too close` |
| `opts.type` | `char extra condition (on top of feature closeness) required before an instance is marked redundant: 'Ftr' (none), 'Ftr&AP' (similar algorithm performance too), 'Ftr&Good' (both instances good on every algorithm), or 'Ftr&AP&Good' (both)` |
## Output Arguments
| Argument | Description |
|---|---|
| `subsetIndex` | `(ninst x 1) logical; true where the instance was found redundant against an earlier, kept instance` |
| `isDissimilar` | `(ninst x 1) logical; false where an instance triggered a feature-space-closeness check against another` |
| `isVISA` | `(ninst x 1) logical; true where two instances were feature-close but did not meet opts.type's extra condition, so neither was marked redundant despite the proximity ("visually important, but not subsetted away")` |
| `unif` | `feature-space uniformity of the retained (non-redundant) subset: 1 minus the coefficient of variation of nearest-neighbour distances (closer to 1 = more evenly spread, closer to 0 = clustered). Previously computed but assigned to an undefined workspace variable (model.data.unif) and discarded; now a real output.` |

## References

- Smith-Miles, K. & Munoz, M.A. (2023). Instance Space Analysis for Algorithm Testing. ACM Computing Surveys, 55(12), Article 255. <https://doi.org/10.1145/3572895>
- Alipour, H., Munoz, M.A. & Smith-Miles, K. (2023). Enhanced instance space analysis for the maximum flow problem. European Journal of Operational Research, 304(2), 411-428. <https://doi.org/10.1016/j.ejor.2022.04.012>

## See Also

[PRELIM](PRELIM.html) | [SIFTED](SIFTED.html) | [InstanceSpace](InstanceSpace.html)
