# buildIS

Thin backward-compatibility wrapper around InstanceSpace.

Preserves the original buildIS(rootdir) calling convention (a plain
function taking a directory and returning the in-memory model struct)
for callers -- notably the MATILDA web platform -- that invoke this
entry point directly. New code should use InstanceSpace directly:

```
obj = InstanceSpace(rootdir);
obj = obj.build();
```

## Input Arguments

| Argument | Description |
|---|---|
| `rootdir` | directory containing metadata.csv and options.json |

## Output Arguments

| Field | Description |
|---|---|
| `model` | the built InstanceSpace object's obj.model (equivalent to calling InstanceSpace(rootdir).build().model directly) |

## References

- Smith-Miles, K. & Munoz, M.A. (2023). Instance Space Analysis for Algorithm Testing. ACM Computing Surveys, 55(12), Article 255. <https://doi.org/10.1145/3572895>

## See Also

[InstanceSpace](InstanceSpace.html)
