# SIFTED

Select the features that best explain algorithm performance

<!-- opts: sifted -->

## Syntax

```
[X,out] = SIFTED(X,Y,Ybin,featlabels,opts)
```

## Description

`[X,out] = SIFTED(X,Y,Ybin,featlabels,opts)` returns the columns of `X` that best explain algorithm performance, and the indices of those columns in `out.selvars`. Selection has two steps:

1. **Correlation filter.** For every algorithm, keep the feature most correlated with its performance. Also keep every feature whose absolute Pearson correlation with the performance of at least one algorithm is at least `opts.rho` and significant at level `opts.pval`.
2. **Clustering and genetic search.** Group the remaining features into `opts.K` clusters by correlation distance (k-means). A genetic algorithm then picks one feature per cluster, scoring each candidate set by how well a *k*-nearest-neighbour classifier predicts `Ybin` in a `PILOT` projection of those features.

Step 2 is skipped when three or fewer features remain, or when no more than `opts.K` remain.

## Examples

### Select features from the reference data

```matlab
rootdir = 'test/data/example/';
if ~isfolder(rootdir), mkdir(rootdir); end
copyfile('test/data/metadata.csv', rootdir);
opts.perf = struct('MaxPerf', false, 'AbsPerf', true, 'epsilon', 0.20);
obj = InstanceSpace(rootdir, opts);
obj = obj.build('stages', {'prelim'});

d = obj.model.data;
siftedOpts = obj.opts.sifted;
siftedOpts.K = 5;
[Xs, out] = SIFTED(d.X, d.Y, d.Ybin, d.featlabels, siftedOpts);
d.featlabels(out.selvars)
```

### Check the choice of K

`out.eva` holds the silhouette value for each number of clusters. When the value for `opts.K` is poor, SIFTED also suggests a better `K`.

```matlab
plot(out.eva.InspectedK, out.eva.CriterionValues, '-o')
xlabel('K'), ylabel('Mean silhouette')
if isfield(out, 'Ksuggested'), disp(out.Ksuggested), end
```

## Input Arguments

### `X` — Feature matrix

*numeric matrix*

Preprocessed features, one row per instance.

### `Y` — Performance matrix

*numeric matrix*

Preprocessed performance, one row per instance.

### `Ybin` — Good-performance labels

*logical matrix*

From `PRELIM`. Used by the genetic search's fitness function.

### `featlabels` — Feature names

*cell array of character vectors*

### `opts` — Selection options

*structure*

Normally `obj.opts.sifted`.

#### `opts.rho` — Minimum correlation

*`0.10` (default) | scalar in [0, 1]*

#### `opts.pval` — Significance level

*`0.05` (default) | scalar in (0, 1)*

#### `opts.K` — Number of clusters

*`10` (default) | positive integer*

Number of features SIFTED selects. When no more than `K` features pass the correlation filter, all of them are kept.

#### `opts.MaxIter` — k-means iterations

*`1000` (default) | positive integer*

#### `opts.Replicates` — k-means replicates

*`100` (default) | positive integer*

#### `opts.seed` — Random seed

*`opts.general.seed` (default) | integer*

Seeds k-means, the cross-validation partition, and the genetic algorithm.

#### `opts.dims` — Projection dimension for the fitness function

*`2` (default) | `3`*

`InstanceSpace` sets this from `opts.pilot.dims`, so features are chosen for the same dimension as the final projection.

## Output Arguments

### `X` — Selected features

*numeric matrix*

The columns `out.selvars` of the input `X`.

### `out` — Selection details

*structure*

#### `out.selvars` — Selected feature indices

Indices into the columns of the input `X`.

#### `out.rho`, `out.p` — Correlations

Correlation of each feature with each algorithm's performance, and the p-values.

#### `out.eva` — Cluster evaluation

`evalclusters` result with the silhouette value for each number of clusters.

#### `out.clust` — Cluster membership

`nfeats`-by-`K` logical matrix.

#### `out.Ksuggested` — Suggested number of clusters

Present only when the silhouette value for `opts.K` is below 0.5.

## Tips

- `opts.sifted.flag = false` skips SIFTED in `InstanceSpace` and keeps every feature.
- SIFTED is the slowest stage. Set `opts.general.parallel = true` to evaluate the genetic algorithm's population in parallel.

## Version History

### v0.9.1 — Seed control

k-means, the partition and the genetic algorithm use `opts.seed`.

### v0.9.0 — Renamed from SIFTED2

`SIFTED2` remains as a deprecated alias.

## References

- Smith-Miles, K. & Muñoz, M.A. (2023). Instance Space Analysis for Algorithm Testing. *ACM Computing Surveys*, 55(12), Article 255. <https://doi.org/10.1145/3572895>

## See Also

`PRELIM` | `PILOT` | `InstanceSpace`
