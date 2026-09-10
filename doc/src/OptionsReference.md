# OptionsReference

Documents every `opts.*` field across all pipeline stages, sourced from `utils/ISAdefaults.m` -- the single source of truth for defaults. Any field you don't set falls back to the default shown here.

## opts.general

| Field        | Default | Description                                                                 |
|--------------|---------|-----------------------------------------------------------------------------|
| `seed`       | 42      | RNG seed                                                                     |
| `verbose`    | true    | Enable verbose output                                                          |
| `parallel`   | false   | Use parallel processing (conservative default)                                |
| `ncores`     | 18      | Number of CPU cores to use for parallel processing                              |

## opts.perf

| Field         | Default | Description                                                                 |
|---------------|---------|-----------------------------------------------------------------------------|
| `MaxPerf`     | false   | Maximize performance instead of minimizing                                    |
| `AbsPerf`     | false   | Use absolute performance threshold instead of relative                      |
| `epsilon`     | 0.05    | Good-performance threshold                                                  |
| `betaThreshold` | 0.55  | Easy-instance fraction (fraction of instances that are easy to solve)       |

## opts.prelim

| Field            | Default | Description                                                                 |
|------------------|---------|-----------------------------------------------------------------------------|
| `iqrMultiplier`  | 5       | Outlier bound = median +/- N*IQR                                             |
| `nanThreshold`   | 0.20    | Maximum fraction of missing values before an instance/feature is dropped     |

## opts.auto

| Field          | Default | Description                                                                 |
|----------------|---------|-----------------------------------------------------------------------------|
| `preproc`      | true    | Run automatic pre-processing                                                   |

## opts.bound

| Field      | Default | Description                                                                 |
|------------|---------|-----------------------------------------------------------------------------|
| `flag`     | true    | Bound outliers                                                                 |

## opts.norm

| Field      | Default | Description                                                                 |
|------------|---------|-----------------------------------------------------------------------------|
| `flag`     | true    | Apply Box-Cox + Z-score normalisation                                         |

## opts.selvars

| Field                | Default | Description                                                                 |
|----------------------|---------|-----------------------------------------------------------------------------|
| `smallscaleflag`     | false   | Enable a small-scale experiment: randomly hold out a fraction of instances (via `cvpartition`) |
| `smallscale`         | 0.30    | Fraction of instances to keep when `smallscaleflag` is true |
| `fileidxflag`        | false   | Enable instance subsetting from an explicit list of instance indices read from a file |
| `fileidx`            | ''      | Path to the file listing instance indices to keep, when `fileidxflag` is true |
| `densityflag`        | false   | Enable density-based instance subsetting via FILTER |
| `mindistance`        | 0.10    | Feature-space distance threshold for FILTER's density-based instance subsetting, when `densityflag` is true |
| `type`               | 'Ftr&Good' | FILTER's redundancy condition: 'Ftr', 'Ftr&AP', 'Ftr&Good', or 'Ftr&AP&Good' |

## opts.sifted

| Field    | Default       | Description                                                                 |
|----------|---------------|-----------------------------------------------------------------------------|
| `flag`   | `true`        | Automated feature selection on/off.                                           |
| `rho`    | `0.10`        | Minimum absolute correlation to keep.                                        |
| `pval`   | `0.05`        | Significance level for correlation.                                           |
| `K`      | `10`          | Number of feature clusters.                                                   |
| `MaxIter`| `1000`        | K-means maximum iterations.                                                   |
| `Replicates` | `100`       | K-means replicates.                                                           |
| `seed`   | `opts.general.seed` | Defaults from the general seed.                                              |

## opts.pilot

| Field        | Default       | Description                                                                 |
|--------------|---------------|-----------------------------------------------------------------------------|
| `analytic`   | `false`       | Use closed-form eigenvector solution when true.                              |
| `ntries`     | `10`          | BFGS multi-start restarts.                                                   |
| `dims`       | `2`           | Projection dimensionality, 2 or 3.                                          |
| `method`     | `'standard'`  | 'standard' BFGS/analytic, or 'pls' Partial Least Squares.                   |
| `alpha`      | `1.0`         | Performance-reconstruction cost weight.                                      |
| `viewGroups` | `{}`          | Per-algorithm-group viewpoints for 3D.                                       |
| `topoWeight` | `0`           | Reserved for future use.                                                     |
| `verbose`    | `opts.general.verbose` | Defaults from general.                                                      |
| `seed`       | `opts.general.seed` | Defaults from general.                                                      |

## opts.cloister

| Field          | Default       | Description                                                                 |
|----------------|---------------|-----------------------------------------------------------------------------|
| `pval`         | `0.05`        | Significance level for feature correlations.                                 |
| `corrThreshold`| `0.70`        | |correlation| above which a corner combination is discarded.                         |
| `maxFeatures`  | `20`          | Feature-count guard before falling back to a plain convex hull.             |

## opts.pythia

| Field            | Default       | Description                                                                 |
|------------------|---------------|-----------------------------------------------------------------------------|
| `flag`           | `true`        | Enable the PYTHIA stage.                                                     |
| `kFold`          | `5`           | Cross-validation folds.                                                     |
| `tuning`         | `'sobol'`     | 'sobol' scrambled Sobol search, 'bayes' MATLAB bayesopt, or pre-supplied params.|
| `nTuningIter`    | `20`          | Sobol/Bayes evaluation budget.                                             |
| `params`         | `[]`          | Pre-calculated hyperparameters, required when tuning='none'.                 |
| `skip`           | `false`       | Bypass training entirely.                                                   |
| `ispolykrnl`     | `false`       | SVM only: polynomial vs Gaussian kernel.                                    |
| `useweights`     | `false`       | Cost-sensitive classification.                                              |
| `ensembleMethod` | `'Bag'`       | fitcensemble method.                                                         |
| `verbose`        | `opts.general.verbose` | Defaults from general.                                                      |
| `seed`           | `opts.general.seed` | Defaults from general.                                                      |
| `classifier`     | `'knn'`       | Registry name: 'knn' default, also 'svm','tree','nb','linear','ensemble'.   |

## opts.trace

| Field          | Default       | Description                                                                 |
|----------------|---------------|-----------------------------------------------------------------------------|
| `method`       | `'trace3'`    | 'trace3' default, or 'legacy' for the pre-refactor DBSCAN+alpha-shape algorithm.|
| `PI`           | `0.6`         | Minimum purity threshold.                                                   |
| `minInstances` | `4`           | Minimum instances for a valid footprint.                                     |
| `minAreaFrac`  | `0.01`        | Footprint must exceed this fraction of space.                                |
| `contra`       | `false` (`true` when `method='legacy'`) | Contradiction removal between best-algorithm footprints.|

## opts.outputs

| Field          | Default       | Description                                                                 |
|----------------|---------------|-----------------------------------------------------------------------------|
| `csv`          | `true`        | Write CSV output files.                                                     |
| `png`          | `true`        | Write PNG figure output files.                                              |
| `fig`          | `true`        | Write a .fig file alongside each 3D footprint PNG.                          |
| `web`          | `false`       | Write colour-scaled CSV data for MATILDA's web tools.                       |

## See Also

[MetadataFormat](MetadataFormat.html) | [MigratingLegacyModel](MigratingLegacyModel.html) | [ISAdefaults](ISAdefaults.html)
