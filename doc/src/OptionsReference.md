# Options Reference

Every field of the options structure, with its default

Options control every stage of the pipeline. Give them to `InstanceSpace` as a structure, or as an `options.json` file in the data folder with the same nesting:

```matlab
opts.perf.MaxPerf = false;
opts.pilot.dims = 3;
obj = InstanceSpace(rootdir, opts);
```

```
{"perf": {"MaxPerf": false}, "pilot": {"dims": 3}}
```

Set only the fields you want to change. `ISAvalidateOpts` checks the fields you set, and `ISAdefaults` fills in the rest with the defaults below. Between staged `build` calls you can change `obj.opts` directly.

## opts.general

Settings for the whole pipeline.

| Field | Default | Description |
|---|---|---|
| `seed` | `42` | Random seed. Every stochastic stage derives its seed from this value, so the same seed and data reproduce a run. |
| `verbose` | `true` | Detailed progress output. Stage start and end messages are always printed. |
| `parallel` | `false` | Use a parallel pool (Parallel Computing Toolbox) for SIFTED, PILOT, PYTHIA and TRACE. |
| `ncores` | `18` | Number of workers of the parallel pool. |

## opts.perf

How performance is judged. Used by `PRELIM`.

| Field | Default | Description |
|---|---|---|
| `MaxPerf` | `false` | `true` if larger performance values are better; `false` for a cost such as error or run time. |
| `AbsPerf` | `false` | `true`: good means better than `epsilon`. `false`: good means within a fraction `epsilon` of the best algorithm on the instance. |
| `epsilon` | `0.05` | Good-performance threshold. In [0, 1] when `AbsPerf` is `false`; any real number otherwise. |
| `betaThreshold` | `0.55` | An instance is beta-easy when more than this fraction of the algorithms are good on it. |

## opts.prelim

Data preparation. Used by `INIT` and `PRELIM`.

| Field | Default | Description |
|---|---|---|
| `iqrMultiplier` | `5` | Features are bounded to `median ± iqrMultiplier*IQR`. |
| `nanThreshold` | `0.20` | A feature with at least this fraction of missing values is removed. |

## opts.auto

Used by `PRELIM`.

| Field | Default | Description |
|---|---|---|
| `preproc` | `true` | Run the automatic outlier bounding and normalisation. `false` leaves the data as given. |

## opts.bound

Used by `PRELIM`.

| Field | Default | Description |
|---|---|---|
| `flag` | `true` | Bound feature outliers. Features with very little variation can have an IQR of zero; remove them or turn bounding off. |

## opts.norm

Used by `PRELIM`.

| Field | Default | Description |
|---|---|---|
| `flag` | `true` | Apply Box-Cox and Z-score transforms, so features and performance are close to normally distributed. Recommended, because PILOT is a linear projection and `CLOISTER` expects centred data. |

## opts.selvars

Which features, algorithms and instances to use. Used by `INIT`, `InstanceSpace` and `FILTER`.

| Field | Default | Description |
|---|---|---|
| `feats` | all | Cell array of feature column names to use, with their `feature_` prefix. |
| `algos` | all | Cell array of algorithm column names to use, with their `algo_` prefix. |
| `smallscaleflag` | `false` | Build from a random fraction of the instances. Useful to try settings on a large dataset. |
| `smallscale` | `0.30` | Fraction of instances kept when `smallscaleflag` is `true`. |
| `fileidxflag` | `false` | Build from the instances listed in the file `fileidx`. |
| `fileidx` | `''` | CSV file with one column of instance indices (row numbers of `metadata.csv`). |
| `densityflag` | `false` | Remove near-duplicate instances with `FILTER`. |
| `mindistance` | `0.10` | `FILTER` distance threshold in feature space. |
| `type` | `'Ftr&Good'` | `FILTER` removal condition: `'Ftr'`, `'Ftr&AP'`, `'Ftr&Good'` or `'Ftr&AP&Good'`. |

## opts.sifted

Feature selection. Used by `SIFTED`.

| Field | Default | Description |
|---|---|---|
| `flag` | `true` | Run SIFTED. `false` keeps every feature. |
| `rho` | `0.10` | Minimum absolute correlation between a feature and an algorithm's performance. |
| `pval` | `0.05` | Significance level of the correlations. |
| `K` | `10` | Number of feature clusters, which is the number of features selected. We recommend 10 or fewer. |
| `MaxIter` | `1000` | Maximum k-means iterations. |
| `Replicates` | `100` | Number of k-means replicates. |
| `seed` | `opts.general.seed` | Random seed of SIFTED. |

## opts.pilot

Projection. Used by `PILOT` and `PILOTviewpoint`.

| Field | Default | Description |
|---|---|---|
| `dims` | `2` | Dimension of the instance space, 2 or 3. The legacy `ISA3D = true` is read as `dims = 3`. |
| `method` | `'standard'` | `'standard'` (analytic or BFGS) or `'pls'` (Partial Least Squares). |
| `analytic` | `false` | Closed-form solution instead of BFGS. Faster, but can be poorly conditioned. |
| `ntries` | `10` | Number of BFGS restarts, and of `PILOTviewpoint` restarts. |
| `alpha` | `1.0` | Weight of performance reconstruction relative to feature reconstruction. Standard method only. |
| `viewGroups` | `{}` | 3D only: cell array of algorithm index vectors; one viewpoint is found per group. Empty means one viewpoint for all algorithms. |
| `topoWeight` | `0` | Reserved; has no effect in this version. |
| `seed` | `opts.general.seed` | Random seed of PILOT and PILOTviewpoint. |
| `verbose` | `opts.general.verbose` | PILOT progress output. |

## opts.cloister

Boundary estimation. Used by `CLOISTER`.

| Field | Default | Description |
|---|---|---|
| `pval` | `0.05` | Significance level of the feature correlations. |
| `corrThreshold` | `0.70` | Correlations stronger than this constrain the boundary. Lower values discard more corners and may fail to give a boundary. |
| `maxFeatures` | `20` | With more features, the convex hull of the instances is used as the boundary instead. |

## opts.pythia

Algorithm selection. Used by `PYTHIA`.

| Field | Default | Description |
|---|---|---|
| `classifier` | `'knn'` | `'knn'`, `'svm'`, `'tree'`, `'nb'`, `'linear'` or `'ensemble'`. See `ISAgetClassifierFcn`. |
| `tuning` | `'sobol'` | Hyperparameter search: `'sobol'`, `'bayes'` or `'none'` (use `params`). |
| `nTuningIter` | `20` | Number of hyperparameter candidates evaluated. |
| `kFold` | `5` | Number of cross-validation folds. |
| `params` | `[]` | Fixed hyperparameters, one row per algorithm and one column per hyperparameter (1 for `'tree'`, `'nb'`, `'linear'`; 2 otherwise). Required with `tuning = 'none'`. |
| `useweights` | `false` | Cost-sensitive training, weighting instances by how far their performance is from the mean. |
| `ispolykrnl` | `false` | SVM only: polynomial kernel instead of Gaussian. |
| `ensembleMethod` | `'Bag'` | `fitcensemble` method for `'ensemble'`. |
| `skip` | `false` | Skip training; `TRACE` then uses the true labels only. |
| `flag` | `true` | Kept for compatibility with older option files; has no effect. |
| `seed` | `opts.general.seed` | Random seed of PYTHIA. |
| `verbose` | `opts.general.verbose` | PYTHIA progress output. |

## opts.trace

Footprints. Used by `TRACE`.

| Field | Default | Description |
|---|---|---|
| `method` | `'trace3'` | `'trace3'`, or `'legacy'` for the earlier DBSCAN method (2D only). |
| `PI` | `0.6` | Minimum purity of a footprint. |
| `minInstances` | `4` | Minimum number of instances in a footprint. |
| `minAreaFrac` | `0.01` | Minimum footprint size as a fraction of the whole space. |
| `contra` | `false` | Legacy method only: remove contradictions between best-algorithm footprints. `true` by default when `method` is `'legacy'`. |

## opts.outputs

Files written when `build` or `explore` completes. See `scriptcsv`, `scriptpng` and `scriptweb`.

| Field | Default | Description |
|---|---|---|
| `csv` | `true` | Write CSV files. |
| `png` | `true` | Write PNG figures. |
| `fig` | `true` | 3D only: also save footprint figures as `.fig` files. |
| `web` | `false` | Write the colour files used by MATILDA. Needs `csv`. |

## See Also

`InstanceSpace` | `ISAdefaults` | `ISAvalidateOpts` | [Metadata File Format](MetadataFormat.html)
