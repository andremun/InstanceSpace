# ISArecallView

Return a 3D figure to its optimised viewpoint

## Syntax

```
ISArecallView(fig)
ISArecallView(fig,algoIdx)
```

## Description

`ISArecallView(fig)` sets the camera of the 3D instance-space figure `fig` to the global viewpoint found by `PILOTviewpoint`. Use it after rotating a figure saved by `scriptpng`.

`ISArecallView(fig,algoIdx)` uses the viewpoint of the group in `opts.pilot.viewGroups` that contains algorithm `algoIdx`, or the global viewpoint if no group contains it.

The viewpoints are stored in the figure's `UserData` by `scriptpng`.

## Examples

### Restore the view of a saved footprint

```matlab
fig = openfig('test/data/example/footprint_KNN.fig');
rotate3d(fig, 'on')        % explore the figure, then:
ISArecallView(fig, 6)      % back to the viewpoint for algorithm 6
```

## Input Arguments

### `fig` — Figure

*figure handle*

A 3D figure saved by `scriptpng`.

### `algoIdx` — Algorithm index

*`[]` (default) | positive integer*

Column of the algorithm in `model.data.algolabels`.

## Version History

### v0.9.0 — Introduced

## See Also

`PILOTviewpoint` | `scriptpng`
