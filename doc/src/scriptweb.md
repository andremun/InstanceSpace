# scriptweb

Write colour-scaled CSV data for MATILDA's web tools to rootdir.

## Syntax

```
scriptweb(container,rootdir)
```

## Description

Only useful when opts.outputs.web=true, i.e. results will be served through matilda.unimelb.edu.au; not needed for local/offline use.

## Input Arguments

| Argument | Description |
|---|---|
| `container` | `struct (model from buildIS/InstanceSpace.build() or a testResults entry from exploreIS/InstanceSpace.explore())` |
| `rootdir` | `string (destination directory; trailing slash required)` |

## Output Arguments

| Field | Description |
|---|---|
| `none` | writes colour-scaled CSV files to rootdir as a side effect (void function) |

## References

- Smith-Miles, K. & Munoz, M.A. (2023). Instance Space Analysis for Algorithm Testing. ACM Computing Surveys, 55(12), Article 255. <https://doi.org/10.1145/3572895>
