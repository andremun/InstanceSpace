# TRACE

Compute algorithm footprints in the instance space

<!-- opts: trace -->

## Syntax

```
out = TRACE(Z,Ybin,Yhat,P,beta,algolabels,opts)
out = TRACE(Z,Ybin,Yhat,P,beta,algolabels,opts,trainedTrace)
```

## Description

`out = TRACE(Z,Ybin,Yhat,P,beta,algolabels,opts)` finds each algorithm's *footprint*: the regions of the instance space where the algorithm is good (`out.good`) and where it is the best algorithm (`out.best`), plus the region of *beta-hard* instances (`out.hard`). For each footprint TRACE reports its area (2D) or volume (3D), density and purity, relative to the whole space.

`out = TRACE(Z,Ybin,Yhat,P,beta,algolabels,opts,trainedTrace)` keeps the footprints of an earlier training call, `trainedTrace`, and recomputes their density and purity with new instances. Footprints are not rebuilt.

The default method, TRACE3, works in 2D and 3D. `opts.method = 'legacy'` selects the earlier DBSCAN-based method, which is 2D only.

## Examples

### Footprints of a trained instance space

```matlab
rootdir = 'test/data/example/';
if ~isfolder(rootdir), mkdir(rootdir); end
copyfile('test/data/metadata.csv', rootdir);
opts.perf = struct('MaxPerf', false, 'AbsPerf', true, 'epsilon', 0.20);
obj = InstanceSpace(rootdir, opts);
obj = obj.build();

disp(obj.model.trace.summary)
obj.plot('footprint', 3)       % good and bad instances, footprint of algorithm 3
```

### Rebuild footprints with a stricter purity

```matlab
m = obj.model;
traceOpts = obj.opts.trace;
traceOpts.PI = 0.75;
out = TRACE(m.pilot.Z, m.data.Ybin, m.pythia.Yhat, m.data.P, m.data.beta, ...
            m.data.algolabels, traceOpts);
```

## Input Arguments

### `Z` — Instance coordinates

*`ninst`-by-2 or `ninst`-by-3 numeric matrix*

### `Ybin` — Good-performance labels

*logical matrix*

### `Yhat` — Predicted labels

*logical matrix | `[]`*

`model.pythia.Yhat`. TRACE3 uses only instances that are both good and predicted good, which removes instances the classifier cannot separate. With `[]` it uses `Ybin` alone and warns (`ISA:TRACE3:noPYTHIA`).

### `P` — Best algorithm per instance

*integer vector*

### `beta` — Beta-easy instances

*logical vector*

`~beta` defines the hard-instance footprint.

### `algolabels` — Algorithm names

*cell array of character vectors*

### `opts` — Footprint options

*structure*

Normally `obj.opts.trace`.

#### `opts.method` — Footprint method

*`'trace3'` (default) | `'legacy'`*

`'legacy'` falls back to `'trace3'` with a warning for a 3D space.

#### `opts.PI` — Purity threshold

*`0.6` (default) | scalar in [0, 1]*

Minimum fraction of good instances inside a footprint.

#### `opts.minInstances` — Minimum instances

*`4` (default) | positive integer*

Fewer instances than this give an empty footprint.

#### `opts.minAreaFrac` — Minimum size

*`0.01` (default) | scalar in [0, 1]*

A footprint smaller than this fraction of the whole space is discarded.

#### `opts.contra` — Remove contradictions

*`false` (default; `true` for `'legacy'`) | logical*

Legacy method only: remove the overlap between best-algorithm footprints.

### `trainedTrace` — Trained footprints

*structure*

`model.trace`. Its presence selects evaluation mode.

## Output Arguments

### `out` — Footprints and summary

*structure*

#### `out.space` — Whole space

The convex hull of all instances: `measure` (area or volume), `measureLabel`, `elements`, `density`, `purity`.

#### `out.good` — Good-performance footprints

Cell array, one footprint per algorithm. Each footprint has `polygon` (an `alphaShape`), `measure`, `measureLabel`, `elements`, `goodElements`, `density` and `purity`.

#### `out.best` — Best-algorithm footprints

Same format as `out.good`.

#### `out.hard` — Beta-hard footprint

#### `out.summary` — Summary table

Cell array: for each algorithm, the area or volume of its good and best footprints (absolute and as a fraction of the space), their density (absolute and normalised) and their purity. Written to `footprint_performance.csv` by `scriptcsv`.

## Algorithms

For each performance vector, TRACE3:

1. Takes the instances that are both good and predicted good by PYTHIA, without duplicates.
2. Builds an `alphaShape` around them at the smallest alpha that encloses them.
3. Accepts the shape when its purity reaches `opts.PI`.
4. Otherwise decreases alpha in 100 steps, which tightens the shape around dense groups of good instances, until the purity reaches `opts.PI` or the shape falls below `opts.minAreaFrac` of the space.

The legacy method clusters the good instances with DBSCAN, builds a polygon per cluster, and optionally removes contradictions between best-algorithm footprints (Muñoz & Smith-Miles, 2017).

## Version History

### v0.9.0 — TRACE3

TRACE3 replaces the DBSCAN method as the default and supports 3D instance spaces. The earlier method remains as `opts.method = 'legacy'`.

## References

- Smith-Miles, K. & Muñoz, M.A. (2023). Instance Space Analysis for Algorithm Testing. *ACM Computing Surveys*, 55(12), Article 255. <https://doi.org/10.1145/3572895>
- Simpson, C., Muñoz, M.A., Kandanaarachchi, S. & Campello, R.J.G.B. (2025). ISA3: A 3-dimensional expansion of Instance Space Analysis. *Machine Learning*, 114, 240. <https://doi.org/10.1007/s10994-025-06871-5>
- Muñoz, M.A. & Smith-Miles, K. (2017). Performance analysis of continuous black-box optimization algorithms via footprints in instance space. *Evolutionary Computation*, 25(4), 529–554. <https://doi.org/10.1162/EVCO_a_00194>

## See Also

`PYTHIA` | `PILOT` | `InstanceSpace` | `scriptpng`
