# ISAdefaults

Fill in missing options with their default values

## Syntax

```
opts = ISAdefaults(opts)
```

## Description

`opts = ISAdefaults(opts)` returns `opts` with every missing field set to its default, so that every stage receives a complete options structure. Fields already present are not changed. The defaults are listed in the [Options Reference](OptionsReference.html).

ISAdefaults also maps a few legacy field names to their current names when the current name is absent: `opts.parallel.flag`/`.ncores` to `opts.general.parallel`/`.ncores`, `opts.pilot.ISA3D` to `opts.pilot.dims`, `opts.cloister.cthres` to `opts.cloister.corrThreshold`, `opts.pythia.cvfolds` to `opts.pythia.kFold`, and `opts.pythia.useknn = false` to `opts.pythia.classifier = 'svm'`.

`InstanceSpace` calls ISAdefaults after `ISAvalidateOpts` when an object is created or loaded.

## Examples

### See every default

```matlab
opts = ISAdefaults(struct());
opts.pilot
```

### Complete a partial options file

```matlab
opts = jsondecode(fileread('options.json'));
opts = ISAdefaults(ISAvalidateOpts(opts));
```

## Input Arguments

### `opts` — Options

*structure*

Any subset of the options, for example as decoded from `options.json`.

## Output Arguments

### `opts` — Complete options

*structure*

## Version History

### v0.9.0 — Introduced

Defaults moved from each stage into one function.

## See Also

`ISAvalidateOpts` | `InstanceSpace` | [Options Reference](OptionsReference.html)
