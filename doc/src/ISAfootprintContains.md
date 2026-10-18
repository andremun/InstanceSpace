# ISAfootprintContains

Test footprint membership with an explicit boundary distance

## Syntax

```
inside = ISAfootprintContains(geometry,queries)
inside = ISAfootprintContains(geometry,queries,tolerance)
```

## Description

`inside = ISAfootprintContains(geometry,queries)` tests whether each query point lies in a footprint, including its boundary. The default distance tolerance is zero.

`inside = ISAfootprintContains(geometry,queries,tolerance)` also includes points within the specified Euclidean distance of a boundary edge in 2D or triangular face in 3D. Hole interiors farther from the boundary remain excluded. The geometry and its area or volume are unchanged.

TRACE uses this helper for footprint metrics and evaluation. Set `opts.trace.boundaryTolerance` to configure the pipeline; TRACE saves the fitted value and reuses it during exploration. The distance is in projected-coordinate units, with no implicit epsilon multiplier, automatic scaling or dependence on other points in the query batch. The membership policy is shared with pyInstanceSpace; floating-point geometry can still differ between implementations.

## Examples

### Include the exact boundary of a square

Run `startup` first when calling utilities directly from a fresh session.

```matlab
startup
region = polyshape([0 0; 1 0; 1 1; 0 1]);
queries = [0.5 0.5; 0.5 0; 0.5 -0.0005; 0.5 -0.002];
inside = ISAfootprintContains(region, queries)
% inside = [true; true; false; false]
```

### Choose an explicit distance tolerance

```matlab
inside = ISAfootprintContains(region, queries, 0.001)
% inside = [true; true; true; false]
```

Choose a positive distance from the uncertainty budget of the projected coordinates. A tolerance is not a guarantee that legacy footprint trimming is invariant under arbitrary coordinate perturbations.

## Input Arguments

### `geometry` — Footprint geometry

*2D `polyshape` | 2D or 3D `alphaShape`*

The fitted region, normally a footprint's `polygon` field. Octave TRACE3 uses the toolkit's compatible alpha-shape class; MATLAB's `polyshape` and legacy TRACE remain MATLAB-only. An empty geometry produces false for every query.

### `queries` — Points to test

*numeric matrix*

One row per point, with two or three columns matching the geometry's dimensionality. Nonfinite query points are outside.

### `tolerance` — Boundary distance

*`0` (default) | finite nonnegative real scalar*

Euclidean distance in projection units. Positive values include points at distance less than or equal to this value from a boundary, including hole boundaries.

## Output Arguments

### `inside` — Membership mask

*logical column vector*

One value per query row. Testing a point alone or alongside other queries uses the same rule.

## Version History

### v0.9.2 — Introduced

Closed-boundary membership with explicit distance tolerance for TRACE construction, legacy contradiction handling and evaluation.

## See Also

`TRACE` | `ISAdefaults` | [Options Reference](OptionsReference.html)
