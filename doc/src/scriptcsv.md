# scriptcsv

Write the results of an instance space to CSV files

## Syntax

```
scriptcsv(container,rootdir)
```

## Description

`scriptcsv(container,rootdir)` writes the coordinates, data, predictions, footprints and summary tables of a trained model or an evaluation result to CSV files in `rootdir`. Coordinate columns are `z_1`, `z_2` and, for a 3D space, `z_3`.

`InstanceSpace.build` and `InstanceSpace.explore` call scriptcsv when `opts.outputs.csv` is `true`.

| File | Contents |
|---|---|
| `coordinates.csv` | instance coordinates `Z` |
| `bounds.csv`, `bounds_prunned.csv` | CLOISTER boundary `Zedge` and `Zecorr` (training only) |
| `feature_raw.csv`, `feature_process.csv` | selected features, before and after preprocessing |
| `algorithm_raw.csv`, `algorithm_process.csv` | performance, before and after preprocessing |
| `algorithm_bin.csv` | good-performance labels `Ybin` |
| `good_algos.csv`, `beta_easy.csv`, `portfolio.csv` | number of good algorithms, beta-easy flag, best algorithm |
| `algorithm_svm.csv`, `portfolio_svm.csv` | PYTHIA predictions and selected algorithm |
| `footprint_<algo>_good.csv`, `footprint_<algo>_best.csv` | footprint boundary vertices; separate regions and holes are separated by a row of `NaN` |
| `footprint_performance.csv` | TRACE summary table |
| `classifier_table.csv` | PYTHIA summary table |
| `projection_matrix.csv` | PILOT projection matrix (training only) |

## Examples

### Write the CSV files of an evaluation to another folder

```matlab
res = obj.getResults(1);
outdir = 'results/test_csv/';
if ~isfolder(outdir), mkdir(outdir); end
scriptcsv(res, outdir);
```

## Input Arguments

### `container` — Model or evaluation result

*structure*

`obj.model` or an element of `obj.testResults`.

### `rootdir` — Output folder

*character vector*

Must exist and end with a file separator.

## Version History

### v0.9.2 — Footprint holes

A footprint region with a hole writes both its outer boundary and the hole's boundary, separated by a `NaN` row.

### v0.9.1 — Multiple footprint regions

Every region of a footprint is written, separated by `NaN` rows, instead of only the first.

## See Also

`scriptpng` | `scriptweb` | `InstanceSpace`
