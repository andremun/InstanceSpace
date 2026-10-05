# ISAmigrateModel

Update a model from an earlier toolbox version to the current layout

## Syntax

```
model = ISAmigrateModel(rootdir)
model = ISAmigrateModel(rootdir,'backupSuffix',suffix)
model = ISAmigrateModel(model)
```

## Description

`model = ISAmigrateModel(rootdir)` migrates `rootdir/model.mat` in place. The original file is first copied to `model_legacy.mat` in the same folder. Use this form once, to convert a saved model.

`model = ISAmigrateModel(rootdir,'backupSuffix',suffix)` names the backup `model<suffix>.mat`. An existing backup is never overwritten; ISAmigrateModel raises `ISA:ISAmigrateModel:backupExists` instead.

`model = ISAmigrateModel(model)` migrates a model structure already in memory and writes no files. `InstanceSpace.load` and `exploreIS` use this form on every load, so a legacy `model.mat` also works without converting it first.

## Examples

### Convert a saved model

```matlab
ISAmigrateModel('/path/to/old/run/');
obj = InstanceSpace.load('/path/to/old/run/');
```

### Migrate in memory

```matlab
m = load('/path/to/old/run/model.mat');
m = ISAmigrateModel(m);
m.pythia.classifiers
```

## Input Arguments

### `rootdir` — Folder containing model.mat

*character vector | string scalar*

### `model` — Model structure

*structure*

A model loaded with `load('model.mat')`.

### `suffix` — Backup file suffix

*`'_legacy'` (default) | character vector*

File form only.

## Output Arguments

### `model` — Migrated model

*structure*

## Algorithms

| Legacy content | Becomes |
|---|---|
| `opts.oracle`, `opts.pbldr`, `opts.sbound`, `opts.footprint` | `opts.pythia`, `opts.pilot`, `opts.cloister`, `opts.trace` |
| `opts.corr.flag`, `opts.corr.threshold`, `opts.clust.flag` | fields of `opts.sifted` |
| `opts.perf.MaxMin` | `opts.perf.MaxPerf` |
| `model.data.bestPerformace` | `model.data.Ybest` |
| PLS model without `pilot.Xmean` | recover the mean from training features and translate saved `cloist.Zedge`/`Zecorr` by `-Xmean*A'`; face indices are unchanged |
| `model.pythia.svm{i}`, `model.pythia.knn{i}` | `model.pythia.classifiers{i}` |
| `model.pythia.boxcosnt`, `model.pythia.kscale` | `model.pythia.param1`, `model.pythia.param2` |
| LIBSVM classifier structures | retrained with `opts.pythia.classifier` (default `'knn'`), because LIBSVM structures have no `predict` method |
| `model.trace` in the DBSCAN/polyshape format | recomputed with TRACE3, using `model.pythia.Yhat` when available, otherwise `model.data.Ybin` |
| no `completedStages` | inferred from the stage fields present |

A `model.pilot.A` without `B` and `C` cannot be repaired and raises a warning.

## GNU Octave

File-based migration recognizes versioned Octave archives and preserves their explicit reconstruction format. MATLAB legacy migration retains its existing file format.

## Version History

### v0.9.2 — PLS boundary alignment and Octave archives

File-based migration now reads versioned Octave model archives. Schema-1 archives remain readable; saving writes schema 2, which keeps ordinary struct fields separate from internal object markers.

Recovering a missing PLS feature mean also translates any saved CLOISTER vertices into the fitted PLS coordinate system. Empty or absent boundaries are retained. Repeated migration leaves already centred geometry unchanged. Models without the training features needed to recover the mean must be rebuilt.

### v0.9.0 — Introduced

## See Also

[Migrating a Legacy Model](MigratingLegacyModel.html) | `InstanceSpace` | `exploreIS`
