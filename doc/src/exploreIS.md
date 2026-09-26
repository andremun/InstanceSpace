# exploreIS

Evaluate a saved instance space on new instances (backward-compatible wrapper)

## Syntax

```
out = exploreIS(rootdir)
```

## Description

`out = exploreIS(rootdir)` loads `rootdir/model.mat`, evaluates `rootdir/metadata_test.csv` with it, writes the output files to `rootdir`, prints `EOF:SUCCESS`, and returns the result. It is equivalent to

```matlab
obj = InstanceSpace.load(rootdir);
obj = obj.explore(rootdir);
out = obj.getResults(1);
```

A `model.mat` from a version before v0.9.0 is migrated in memory by `ISAmigrateModel`. For new code, use `InstanceSpace`.

## Examples

### Evaluate test instances after buildIS

```matlab
copyfile('test/data/metadata_test.csv', rootdir);
model = buildIS(rootdir);
out = exploreIS(rootdir);
out.pythia.summary
```

## Input Arguments

### `rootdir` — Data folder

*character vector*

Folder that contains `model.mat` and `metadata_test.csv`.

## Output Arguments

### `out` — Evaluation result

*structure*

Same fields as a trained model, computed for the test instances.

## Version History

### v0.9.0 — Wrapper around InstanceSpace

`exploreIS` became a wrapper around `InstanceSpace`. Its inputs, outputs and files are unchanged.

## See Also

`buildIS` | `InstanceSpace` | `ISAmigrateModel`
