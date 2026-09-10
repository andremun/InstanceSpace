# ISArecallView

Snap an open 3D instance-space figure to its stored PILOTviewpoint camera angle for a given algorithm group.

## Syntax

```
ISArecallView(fig, groupIdx)
ISArecallView(fig)          % use the default/global viewpoint
```

## Description

Useful after manually rotating a 3D footprint .fig while exploring it interactively, to return to the optimised viewpoint without having to recompute or look it up by hand.

## Input Arguments

| Argument | Description |
|---|---|
| `fig` | handle to a figure (produced by scriptpng.m or a .fig file) containing viewpoint data in UserData. |
| `groupIdx` | algorithm column index (from opts.pilot.viewGroups) to apply a specific viewpoint, or [] for the default/global viewpoint. |

## Output Arguments

| Field | Description |
|---|---|
| `None` | rotates the camera to the stored viewpoint as a side effect. |

## References

- Smith-Miles, K. & Munoz, M.A. (2023). Instance Space Analysis for Algorithm Testing. *ACM Computing Surveys*, 55(12), Article 255. <https://doi.org/10.1145/3572895>
- Simpson, C., Munoz, M.A., Kandanaarachchi, S. & Campello, R.J.G.B. (2025). ISA3: A 3-dimensional expansion of Instance Space Analysis. *Machine Learning*, 114, 240. <https://doi.org/10.1007/s10994-025-06871-5>
- Munoz, M.A., Villanova, L., Baatar, D. & Smith-Miles, K. (2018). Instance spaces for machine learning classification. *Machine Learning*, 107(1), 109-147. <https://doi.org/10.1007/s10994-017-5629-5>

## See Also

[scriptpng](scriptpng.html) | [PILOTviewpoint](PILOTviewpoint.html)
