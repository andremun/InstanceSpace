# scriptpng

Save figures of an instance space as PNG files

## Syntax

```
scriptpng(container,rootdir)
```

## Description

`scriptpng(container,rootdir)` draws the standard views of a trained model or an evaluation result and saves them in `rootdir`: each feature, each algorithm's performance, good/bad labels, predictions, footprints, the algorithm portfolio, the number of good algorithms, the beta score, instance sources, and the CLOISTER boundary.

For a 3D space the figures use the viewpoints found by `PILOTviewpoint`. Each footprint is also saved as a MATLAB figure (`footprint_<algo>.fig`, `footprint_portfolio.fig`) for interactive rotation, unless `opts.outputs.fig` is `false`. `ISArecallView` returns a rotated figure to its stored viewpoint.

`InstanceSpace.build` and `InstanceSpace.explore` call scriptpng when `opts.outputs.png` is `true`.

| File | Figure |
|---|---|
| `distribution_feature_<feat>.png` | feature value over the space |
| `distribution_performance_global_normalized_<algo>.png`, `..._individual_normalized_<algo>.png` | performance, scaled over all algorithms or per algorithm |
| `binary_performance_<algo>.png`, `binary_classifier_<algo>.png` | good/bad labels, and PYTHIA's predictions |
| `footprint_<algo>.png`, `footprint_portfolio.png` | footprints |
| `distribution_portfolio.png`, `distribution_svm_portfolio.png` | best and selected algorithm |
| `distribution_number_good_algos.png`, `distribution_beta_score.png` | number of good algorithms, beta-easy instances |
| `distribution_sources.png` | instance sources, if the metadata has a `source` column |
| `distribution_boundary.png` | CLOISTER boundary (training only) |

## Examples

### Save figures at a larger font size

```matlab
outdir = 'results/figures/';
if ~isfolder(outdir), mkdir(outdir); end
set(groot, 'defaultAxesFontSize', 14)
scriptpng(obj.model, outdir);
```

## Input Arguments

### `container` — Model or evaluation result

*structure*

`obj.model` or an element of `obj.testResults`.

### `rootdir` — Output folder

*character vector*

Must exist and end with a file separator.

## Version History

### v0.9.2 — 3D boundary figure

`distribution_boundary.png` is also drawn for a 3D space.

### v0.9.1 — Boundary figure

Adds `distribution_boundary.png`.

### v0.9.0 — 3D figures

3D spaces are drawn from their optimised viewpoints and saved as `.fig` files too.

## See Also

`scriptcsv` | `ISArecallView` | `PILOTviewpoint` | `InstanceSpace`
