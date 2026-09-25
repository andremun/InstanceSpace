# INIT

Read and filter instance metadata

## Syntax

```
[data,extra] = INIT(rootdir,opts)
[data,extra] = INIT(rootdir,opts,trainedModel)
```

## Description

`[data,extra] = INIT(rootdir,opts)` reads `rootdir/metadata.csv`, keeps the features and algorithms listed in `opts.selvars.feats` and `opts.selvars.algos` (all of them when these fields are absent), and returns the matrices the other stages need.

`[data,extra] = INIT(rootdir,opts,trainedModel)` reads `rootdir/metadata_test.csv` for evaluation. It checks that the test file has the trained model's features in the same order, and aligns the algorithm columns with the trained model: known algorithms go to their trained position, new algorithms are appended, and trained algorithms missing from the test file get a column of `NaN`.

## Examples

### Read the reference metadata

```matlab
opts = ISAdefaults(struct());
data = INIT('test/data/', opts);
size(data.X)          % instances x features
data.algolabels       % algorithm names, without the algo_ prefix
```

### Read only some algorithms

```matlab
opts.selvars.algos = {'algo_NB', 'algo_KNN', 'algo_RandF'};
data = INIT('test/data/', opts);
```

## Input Arguments

### `rootdir` — Data folder

*character vector*

Folder that contains `metadata.csv` (training) or `metadata_test.csv` (evaluation). Must end with a file separator.

### `opts` — Toolbox options

*structure*

The complete options structure. INIT reads `opts.selvars.feats` and `opts.selvars.algos`, cell arrays of column names including their `feature_`/`algo_` prefix.

### `trainedModel` — Trained model

*structure*

The `model` property of a built `InstanceSpace` object. Its presence selects evaluation mode.

## Output Arguments

### `data` — Metadata

*structure*

#### `data.instlabels` — Instance names

From the `Instances` column.

#### `data.X`, `data.Xraw` — Features

Numeric matrix; both are the unprocessed values at this point.

#### `data.Y`, `data.Yraw` — Performance

Numeric matrix; both are the unprocessed values at this point.

#### `data.algolabels` — Algorithm names

Without the `algo_` prefix.

#### `data.featlabels` — Feature names

Training mode only, without the `feature_` prefix.

#### `data.S` — Instance sources

Categorical vector, present only when the file has a `source` column.

### `extra` — Evaluation bookkeeping

*structure*

Empty in training mode. In evaluation mode: `featlabelsAll`, `modelalgos` (number of trained algorithms), `newalgos` (algorithms new in the test file), `nalgos` and `ninst`.

## Tips

- See [Metadata File Format](MetadataFormat.html) for the file layout.
- Evaluation fails with `ISA:InstanceSpace:featureOrderMismatch` when `metadata_test.csv` lists its features in a different order from the training file.

## Version History

### v0.9.1 — Introduced

Replaces two separate readers inside `InstanceSpace` (build and explore) with one function.

## References

- Smith-Miles, K. & Muñoz, M.A. (2023). Instance Space Analysis for Algorithm Testing. *ACM Computing Surveys*, 55(12), Article 255. <https://doi.org/10.1145/3572895>

## See Also

`PRELIM` | `InstanceSpace` | [Metadata File Format](MetadataFormat.html)
