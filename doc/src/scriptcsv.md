# scriptcsv

Write a model or explore() result to CSV files in rootdir.

## Syntax

```
scriptcsv(container,rootdir)
```

## Description

Writes projected coordinates, feature/performance tables, algorithm selections, and footprint boundary points, sized for 2D or 3D projections according to size(container.pilot.Z,2).

## Input Arguments

| Argument | Description |
|---|---|
| `container` | model struct from buildIS/InstanceSpace.build() or a testResults entry from exploreIS/InstanceSpace.explore(). |
| `rootdir` | destination directory (trailing slash required). |

## Output Arguments

| Field | Description |
|---|---|
| `none` | writes CSV files to rootdir as a side effect. |

## References

- Smith-Miles, K. & Munoz, M.A. (2023). Instance Space Analysis for Algorithm Testing. *ACM Computing Surveys*, 55(12), Article 255. <https://doi.org/10.1145/3572895>

## See Also

[scriptpng](scriptpng.html) | [scriptweb](scriptweb.html) | [InstanceSpace](InstanceSpace.html)
