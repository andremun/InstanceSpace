# buildIS

Build an instance space from a folder (backward-compatible wrapper)

## Syntax

```
model = buildIS(rootdir)
```

## Description

`model = buildIS(rootdir)` builds an instance space from `rootdir/metadata.csv` and `rootdir/options.json`, writes `model.mat` and the output files to `rootdir`, prints `EOF:SUCCESS`, and returns the trained model. It is equivalent to

```matlab
obj = InstanceSpace(rootdir);
obj = obj.build();
model = obj.model;
```

`buildIS` keeps the calling convention of versions before v0.9.0, which the MATILDA web platform and other scripts use. For new code, use `InstanceSpace`.

## Examples

### Build from a folder

```matlab
rootdir = 'test/data/example/';
if ~isfolder(rootdir), mkdir(rootdir); end
copyfile('test/data/metadata.csv', rootdir);
opts.perf = struct('MaxPerf', false, 'AbsPerf', true, 'epsilon', 0.20);
fid = fopen(fullfile(rootdir, 'options.json'), 'w');
fprintf(fid, '%s', jsonencode(opts));
fclose(fid);

model = buildIS(rootdir);
```

## Input Arguments

### `rootdir` — Data folder

*character vector*

Folder that contains `metadata.csv` and, optionally, `options.json`.

## Output Arguments

### `model` — Trained model

*structure*

The `model` property of the built `InstanceSpace` object.

## Version History

### v0.9.0 — Wrapper around InstanceSpace

`buildIS` became a wrapper around `InstanceSpace`. Its inputs, outputs and files are unchanged.

## See Also

`exploreIS` | `InstanceSpace`
