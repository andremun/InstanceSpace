# PILOTviewpoint
Find the optimal 2D camera viewpoint(s) of a 3D PILOT projection (Equation 2 of Simpson et al., 2025).
## Syntax
```
out = PILOTviewpoint(Z, Y, opts)
```
## Description
For each group, jointly fits A (2xn, the view) and C (nalgos x 2, the performance reconstruction) via BFGS, minimising `||Y_group - (C*A*Z')'||_F^2 + LAMBDA*|dot(v1,v2)|` where `v1=A(1,:)`, `v2=A(2,:)`, following the same multi-start / topological-preservation trial-selection scheme as PILOT.m's numerical branch (Hd = pdist(Z), best trial = highest corr(Hd, pdist(Z*A'))). LAMBDA=0.2 is the paper-calibrated orthogonality penalty weight (not user-exposed). v1 and v2 are rescaled to unit magnitude once per trial (both for the topological-preservation scoring and the stored solution), replacing the original PILOTANGLE.m reference's redundant and axis-wrong sum(A.^2') row-normalisation applied twice.
## Input Arguments
| Argument | Description |
|---|---|
| `Z` | (ninst x n) 3D PILOT projection, n==3 (opts.pilot.dims == 3) |
| `Y` | (ninst x nalgos) performance matrix |
| `opts` | struct with fields: |
| `opts.viewGroups` | cell array of algorithm index vectors, one viewpoint computed per group (default {} -> one global viewpoint over all algorithms) |
| `opts.ntries` | BFGS multi-start restarts (default 10, same convention as opts.pilot.ntries) |
| `opts.X0` | optional (2n+2*numel(group)) x ntries user-supplied starting points, reused for every group whose size matches |
## Output Arguments
| Argument | Description |
|---|---|
| `out` | struct with fields: |
| `out.groups` | the resolved cell array of algorithm groups |
| `out.A` | (ngroups x 1) cell array of fitted 2xn view matrices [row1=v1; row2=v2] flattening Z onto the viewing plane |
| `out.azimuth` | (ngroups x 1) azimuth angle (radians) of the viewing direction cross(v1,v2), for view(az,el) |
| `out.elevation` | (ngroups x 1) elevation angle (radians) |

## References

- Smith-Miles, K. & Munoz, M.A. (2023). Instance Space Analysis for Algorithm Testing. ACM Computing Surveys, 55(12), Article 255. <https://doi.org/10.1145/3572895>
- Simpson, C., Munoz, M.A., Kandanaarachchi, S. & Campello, R.J.G.B. (2025). ISA3: A 3-dimensional expansion of Instance Space Analysis. Machine Learning, 114, 240. <https://doi.org/10.1007/s10994-025-06871-5>
