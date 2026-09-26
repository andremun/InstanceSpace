# PILOTviewpoint

Find the best viewing angles of a 3D instance space

<!-- opts: pilot -->

## Syntax

```
out = PILOTviewpoint(Z,Y,opts)
```

## Description

`out = PILOTviewpoint(Z,Y,opts)` finds, for each group of algorithms in `opts.viewGroups`, the 2D view of the 3D projection `Z` from which the performance of those algorithms is best explained. Each view is a pair of orthogonal directions `v1`, `v2`; the returned azimuth and elevation point the MATLAB camera along `cross(v1,v2)` (Simpson et al., 2025, Equation 2).

`InstanceSpace` calls PILOTviewpoint after PILOT when `opts.pilot.dims` is 3 and stores the result in `model.pilot.viewpoint`. `scriptpng` then renders 3D figures from these angles.

## Examples

### View a 3D instance space from its best angle

```matlab
rootdir = 'test/data/example/';
if ~isfolder(rootdir), mkdir(rootdir); end
copyfile('test/data/metadata.csv', rootdir);
opts.perf = struct('MaxPerf', false, 'AbsPerf', true, 'epsilon', 0.20);
opts.pilot.dims = 3;
obj = InstanceSpace(rootdir, opts);
obj = obj.build('stages', {'prelim', 'sifted', 'pilot'});

Z = obj.model.pilot.Z;
vp = PILOTviewpoint(Z, obj.model.data.Y, obj.opts.pilot);
scatter3(Z(:,1), Z(:,2), Z(:,3), 12, 'filled')
view(rad2deg(vp.azimuth(1)), rad2deg(vp.elevation(1)))
```

### One viewpoint per group of algorithms

```matlab
pilotOpts = obj.opts.pilot;
pilotOpts.viewGroups = {[1 2 3], [4 5], 6:10};
vp = PILOTviewpoint(Z, obj.model.data.Y, pilotOpts);
rad2deg([vp.azimuth vp.elevation])
```

## Input Arguments

### `Z` — 3D instance coordinates

*`ninst`-by-3 numeric matrix*

`model.pilot.Z` from a 3D build.

### `Y` — Performance matrix

*numeric matrix*

Same number of rows as `Z`.

### `opts` — Viewpoint options

*structure*

Normally `obj.opts.pilot`.

#### `opts.viewGroups` — Algorithm groups

*`{}` (default) | cell array of index vectors*

One viewpoint is computed per group, from the columns of `Y` in that group. Empty means one viewpoint for all algorithms.

#### `opts.ntries` — Number of BFGS restarts

*`10` (default) | positive integer*

#### `opts.seed` — Random seed

*`opts.general.seed` (default) | integer*

#### `opts.X0` — Starting points

*numeric matrix*

Optional, `2*3 + 2*numel(group)` rows by one column per trial. Used for each group of matching size.

## Output Arguments

### `out` — Viewpoints

*structure*

#### `out.groups` — Algorithm groups

The groups used, one cell per viewpoint.

#### `out.A` — View matrices

Cell array; each cell is a 2-by-3 matrix `[v1; v2]` that flattens `Z` onto the viewing plane.

#### `out.azimuth`, `out.elevation` — Camera angles

In radians, one per group. Use `view(rad2deg(az), rad2deg(el))`.

## Algorithms

For each group, PILOTviewpoint jointly fits the view `A` (2-by-3) and a performance reconstruction `C` by BFGS, minimising

```
||Y_group - (C*A*Z')'||^2 + 0.2*|dot(v1, v2)|
```

The second term keeps the two view directions orthogonal. Like PILOT, it keeps the trial whose view best preserves the pairwise distances of `Z`.

## Version History

### v0.9.0 — Introduced

Part of the 3D instance space (ISA3) support.

## References

- Simpson, C., Muñoz, M.A., Kandanaarachchi, S. & Campello, R.J.G.B. (2025). ISA3: A 3-dimensional expansion of Instance Space Analysis. *Machine Learning*, 114, 240. <https://doi.org/10.1007/s10994-025-06871-5>

## See Also

`PILOT` | `ISArecallView` | `scriptpng`
