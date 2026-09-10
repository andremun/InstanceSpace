# scriptpng

Write a model or explore() result to PNG figures in rootdir.

## Syntax

```
scriptpng(container, rootdir)
```

## Description

Produces per-feature and per-algorithm distribution plots, portfolio selection and footprint plots, using scriptfcn.m's drawing helpers. Renders in 3D and applies the optimised camera viewpoint (container.pilot.viewpoint, see PILOTviewpoint.m) when the projection is 3D and one was computed. For 3D projections, also writes a .fig file alongside each footprint PNG (footprint_<algo>.fig, footprint_portfolio.fig) for interactive rotation in MATLAB, unless container.opts.outputs.fig is false. Every figure carries the viewpoint struct in its UserData so ISArecallView can snap a reopened .fig back to its optimised camera angle later.

## Input Arguments

| Argument | Description |
|---|---|
| `container` | `struct` (model from `buildIS/InstanceSpace.build()` or a `testResults` entry from `exploreIS/InstanceSpace.explore()`). |
| `rootdir` | `string` (destination directory; trailing slash required). |

## Output Arguments

| Field | Description |
|---|---|
| `none` | writes PNG (and, for 3D projections, .fig) files to rootdir as a side effect (void function). |

## References

- Smith-Miles, K. & Munoz, M.A. (2023). Instance Space Analysis for Algorithm Testing. ACM Computing Surveys, 55(12), Article 255. <https://doi.org/10.1145/3572895>
- Simpson, C., Munoz, M.A., Kandanaarachchi, S. & Campello, R.J.G.B. (2025). ISA3: A 3-dimensional expansion of Instance Space Analysis. Machine Learning, 114, 240. <https://doi.org/10.1007/s10994-025-06871-5>
