# Migrating a Legacy Model

Use a model.mat created by a toolbox version before v0.9.0

Version 0.9.0 renamed several option groups and model fields and replaced the LIBSVM classifiers and the footprint algorithm. A `model.mat` saved by an earlier version therefore needs migrating before the current code can use it. In most cases this is automatic.

## Automatic Migration on Load

`InstanceSpace.load` and `exploreIS` migrate a legacy model in memory each time they load it. Nothing is written to disk:

```matlab
obj = InstanceSpace.load('/path/to/old/run/');   % migrated in memory
obj = obj.explore('/path/to/new/instances/');
```

## Converting the File Once

To avoid repeating the migration, including any classifier retraining, on every load, convert the file:

```matlab
ISAmigrateModel('/path/to/old/run/');
```

The original is kept as `model_legacy.mat`. Choose another backup name with `ISAmigrateModel(rootdir,'backupSuffix','_v08')`.

## What Changes

- **Options** are renamed to the current groups: `opts.oracle` → `opts.pythia`, `opts.pbldr` → `opts.pilot`, `opts.sbound` → `opts.cloister`, `opts.footprint` → `opts.trace`, and `opts.perf.MaxMin` → `opts.perf.MaxPerf`. The feature-selection flags of `opts.corr` and `opts.clust` move to `opts.sifted`.
- **Classifiers.** LIBSVM models cannot be evaluated by the current code. They are retrained with the classifier in `opts.pythia.classifier` (`'knn'` if unset) on the training data stored in the model. Predictions of the retrained model can differ slightly from those of the original.
- **Footprints** in the old DBSCAN/polyshape format are recomputed with TRACE3. Their areas, densities and purities can differ from the values reported by the old version.
- **Field names** such as `data.bestPerformace` are corrected.

The complete table is on the `ISAmigrateModel` reference page.

## When Migration Cannot Complete

- A projection without its reconstruction matrices (`model.pilot.B`, `model.pilot.C`) is kept as is, with a warning.
- Retraining LIBSVM classifiers and recomputing footprints need the training data stored in the model (`model.pilot.Z` and `model.data`). When it is missing, the fields are renamed only, with a warning.
- If a LIBSVM model has to be evaluated as-is, `PYTHIA` raises `ISA:PYTHIA:noLibsvm` unless the LIBSVM MEX files are on the path. The toolbox does not include them; get them from the [LIBSVM project](https://www.csie.ntu.edu.tw/~cjlin/libsvm/). Retraining is recommended instead.

## See Also

`ISAmigrateModel` | `InstanceSpace` | [What's New](WhatsNew.html)
