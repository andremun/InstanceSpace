# MigratingLegacyModel

## Introduction

This guide is designed as a companion to the `ISAmigrateModel` reference page, providing a user-friendly explanation of when and why to migrate a legacy model. Specifically, this guide is relevant if you are using a `model.mat` file that was produced by a pre-v0.9.0 version of the toolkit (before the introduction of the InstanceSpace class). The goal of this migration process is to make these legacy models compatible with the current toolkit's `exploreIS` or `InstanceSpace.explore()` functions, or with output scripts and pipelines.

## When to Migrate

You should consider migrating your legacy model when you are using a `model.mat` file that was created by a version of the toolkit prior to v0.9.0. This migration is necessary to ensure compatibility with the current version of the toolkit, which uses a different class-based structure.

## Migration Methods

There are two primary methods for migrating your legacy model:

### File-based (Recommended)

The recommended approach is to use the `ISAmigrateModel(rootdir)` function. This function reads the `model.mat` file located in `rootdir/model.mat`, backs up the original file to `model_legacy.mat`, and writes the migrated model back to `model.mat`.

**Example:**
```matlab
rootdir = '/path/to/model/directory';
ISAmigrateModel(rootdir);
```

### In-memory

If your legacy model is already loaded into the MATLAB workspace as a struct, you can use the `ISAmigrateModel(model)` function. This method performs the migration without any file I/O operations.

**Example:**
```matlab
model = load('path/to/model.mat');
migratedModel = ISAmigrateModel(model);
```

## High-Level Summary of Changes

The migration process involves several changes to the legacy model's structure. A full field-by-field table is available on the `ISAmigrateModel` reference page. Below is a high-level summary of the changes:

- **Field renames:**
  - `opts` struct updates:
    - `opts.oracle`/`opts.pbldr`/`opts.sbound`/`opts.footprint` -> `opts.pythia`/`opts.pilot`/`opts.cloister`/`opts.trace`
  - `opts` merges:
    - `opts.corr.flag/.threshold` and `opts.clust.flag` merged into `opts.sifted`
  - `opts` updates:
    - `opts.perf.MaxMin` -> `opts.perf.MaxPerf`
  - `model.data` updates:
    - `model.data.bestPerformace` (typo) -> `model.data.Ybest`
  - `model.pilot` updates:
    - `model.pilot.A` without `B/C` -> warning only (not expected; not auto-fixable)
  - `model.pythia` updates:
    - `model.pythia.svm{i}` / `.knn{i}` -> `model.pythia.classifiers{i}`
    - `model.pythia.boxcosnt` / `.kscale` -> `model.pythia.param1` / `.param2`
  - `model.pythia` LIBSVM struct:
    - Retrained via the current classifier registry (default 'knn')
  - `model.trace` updates:
    - Recomputed fresh via TRACE3, using `model.pythia.Yhat` when available (else `model.data.Ybin`)
  - `model.completedStages` updates:
    - Inferred from which sub-structs are present

After migration, the model can be passed to `PYTHIA` eval mode and `scriptcsv` -- i.e., it becomes usable with the current toolkit's `explore()`/output pipeline.

## See Also

For the full field-by-field migration table, see [ISAmigrateModel](ISAmigrateModel.html).
