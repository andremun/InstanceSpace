# ISAsubsetData

Keep a subset of the instances in a data structure

## Syntax

```
data = ISAsubsetData(data,subsetIndex)
data = ISAsubsetData(data,subsetIndex,featIdx)
```

## Description

`data = ISAsubsetData(data,subsetIndex)` keeps the rows `subsetIndex` of every per-instance field of `data`: `X`, `Y`, `Xraw`, `Yraw`, `Ybin`, `beta`, `numGoodAlgos`, `Ybest`, `P`, `instlabels` and, if present, `S`.

`data = ISAsubsetData(data,subsetIndex,featIdx)` also keeps only the feature columns `featIdx` of `data.X` and `data.featlabels`. `InstanceSpace` uses this form after `SIFTED` when density filtering is on.

## Examples

### Keep the beta-easy instances

```matlab
obj = obj.build('stages', {'prelim'});
easy = ISAsubsetData(obj.model.data, obj.model.data.beta);
size(easy.X, 1)
```

## Input Arguments

### `data` — Instance data

*structure*

`model.data`, with the fields listed above.

### `subsetIndex` — Instances to keep

*logical vector | vector of indices*

### `featIdx` — Features to keep

*vector of indices*

## Output Arguments

### `data` — Subset

*structure*

## Version History

### v0.9.0 — Introduced

## See Also

`FILTER` | `InstanceSpace`
