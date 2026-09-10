# SIFTED
Automated feature selection for instance-space analysis.
## Syntax
```
[X, out] = SIFTED(X, Y, Ybin, featlabels, opts)
```
## Description
Selects features in two stages:
1. Correlation filter  -- keeps features whose Pearson correlation with any algorithm's performance exceeds opts.rho and is statistically significant at the opts.pval level.
2. Correlation clustering + GA  -- groups remaining features by correlation distance into opts.K clusters, then uses a genetic algorithm with a KNN fitness function to pick one representative feature per cluster.
## Input Arguments
| Argument | Description |
|---|---|
| `X` | (ninst x nfeats) feature matrix |
| `Y` | (ninst x nalgos) performance matrix |
| `Ybin` | (ninst x nalgos) logical good-performance matrix |
| `featlabels` | (1 x nfeats) cell array of feature name strings |
| `opts.rho` | double minimum absolute correlation to keep (0.10) |
| `opts.pval` | double significance level for correlation (0.05) |
| `opts.K` | int number of feature clusters (10) |
| `opts.MaxIter` | int k-means maximum iterations (1000) |
| `opts.Replicates` | int k-means replicates (100) |
## Output Arguments
| Argument | Description |
|---|---|
| `X` | (ninst x nselected) reduced feature matrix |
| `out.selvars` | indices of selected features (into original X columns) |
| `out.rho`, `out.p` | full correlation and p-value matrices |
| `out.eva` | evalclusters result object |
| `out.clust` | (nfeats x K) cluster membership matrix |
| `out.Ksuggested` | (optional) better K value if current K gives poor silhouette |
## References

- Smith-Miles, K. & Munoz, M.A. (2023). Instance Space Analysis for Algorithm Testing. ACM Computing Surveys, 55(12), Article 255. <https://doi.org/10.1145/3572895>

## See Also

[PRELIM](PRELIM.html) | [PILOT](PILOT.html) | [InstanceSpace](InstanceSpace.html)
