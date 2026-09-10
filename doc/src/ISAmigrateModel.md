# ISAmigrateModel

Migrate a legacy ISA model to the current field layout.

Two calling conventions, dispatched on the type of the first argument:

1) File-based (the primary/recommended form): pass a rootdir
containing model.mat. The original file is backed up alongside it
(default name: model_legacy.mat) and the migrated model is written
back to model.mat in the same directory.

```
ISAmigrateModel(rootdir)
ISAmigrateModel(rootdir, 'backupSuffix', '_v1')   % -> model_v1.mat
```

2) In-memory: pass an already-loaded model struct and use the returned,
migrated struct directly; no file I/O occurs. This form exists
because exploreIS.m already holds the loaded model in memory (it
calls load(modelfile) itself to migrate-then-fill-defaults in one
pass) and has no reason to write the migrated model back to disk
and re-read it on every explore() call — the file-based form is for
one-time offline migration of a model.mat produced by an older
toolkit version, not for use on every read.

```
model = ISAmigrateModel(model);
```

Both paths apply the complete legacy migration table below:

- opts struct renames  
  opts.oracle/opts.pbldr/opts.sbound/opts.footprint  
  -> opts.pythia/opts.pilot/opts.cloister/opts.trace
- opts merges  
  opts.corr.flag/.threshold and opts.clust.flag  
  merged into opts.sifted
- opts.perf.MaxMin -> opts.perf.MaxPerf
- model.data.bestPerformace (typo) -> model.data.Ybest
- model.pilot.A without B/C -> warning only (not expected; not auto-fixable)
- model.pythia.svm{i} / .knn{i} -> model.pythia.classifiers{i}
- model.pythia.boxcosnt / .kscale -> model.pythia.param1 / .param2
- LIBSVM struct in model.pythia -> retrained via the current classifier
  registry (opts.pythia.classifier, default 'knn'); a LIBSVM struct has no
  predict() method, so unlike the plain field renames above it cannot simply
  be relabelled
- model.trace in the pre-refactor DBSCAN+polyshape triangulation format  
  -> recomputed fresh via TRACE3, using model.pythia.Yhat when available (else
  model.data.Ybin)
- missing model.completedStages -> inferred from which sub-structs are
  present (isfield only, same approach as InstanceSpace.load())

After migration the model can be passed to PYTHIA eval mode and scriptcsv.

Examples:
```
ISAmigrateModel('/path/to/rootdir/');            % migrates model.mat on disk
m = load('model.mat'); m = ISAmigrateModel(m);   % migrates an in-memory struct

## Input Arguments

| Argument | Description |
|---|---|
| `input` | either a rootdir (char/string) containing model.mat (file-based form), or an already-loaded model struct (in-memory form); dispatched on this argument's type. |
| `varargin` | name-value pairs, optional; file-based form only. backupSuffix (char, default '_legacy') - suffix appended to the backup filename; must produce a filename different from model.mat, and refuses to overwrite an existing backup at that name. |

## Output Arguments

| Argument | Description |
|---|---|
| `modelOut` | the migrated model struct. File-based form also writes it back to model.mat in rootdir (after backing up the original) as a side effect; in-memory form performs no file I/O -- only the returned struct is migrated. |

## References

- Smith-Miles, K. & Munoz, M.A. (2023). Instance Space Analysis for Algorithm Testing. ACM Computing Surveys, 55(12), Article 255. <https://doi.org/10.1145/3572895>

## See Also

[InstanceSpace](InstanceSpace.html) | [MigratingLegacyModel](MigratingLegacyModel.html)
