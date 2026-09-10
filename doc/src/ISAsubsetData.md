# ISAsubsetData

Subset rows (and optionally feature columns) of a data struct.

## Syntax

```
data = ISAsubsetData(data, subsetIndex)
data = ISAsubsetData(data, subsetIndex, featIdx)
```

## Description

subsets all row-indexed fields. also selects feature columns featIdx from data.X (used in the post-SIFTED density path).

## Input Arguments

| Argument | Description |
|---|---|
| `data` | struct containing fields X, Y, Xraw, Yraw, Ybin, beta, numGoodAlgos, Ybest, P, instlabels, and optionally S. |
| `subsetIndex` | index or logical vector for row selection. |
| `featIdx` | optional column indices for data.X and data.featlabels. |

## Output Arguments

| Field | Description |
|---|---|
| `data` | the subsetted struct with corresponding fields and columns updated. |

## References

- Smith-Miles, K. & Munoz, M.A. (2023). Instance Space Analysis for Algorithm Testing. ACM Computing Surveys, 55(12), Article 255. <https://doi.org/10.1145/3572895>

## See Also

[InstanceSpace](InstanceSpace.html)
