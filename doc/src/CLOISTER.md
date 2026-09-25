# CLOISTER

Estimate the boundary of the instance space

<!-- opts: cloister -->

## Syntax

```
out = CLOISTER(X,A,opts)
```

## Description

`out = CLOISTER(X,A,opts)` estimates the region of the instance space where instances could exist, using correlations between features (Correlated Limits of the Instance Space's Theoretical or Experimental Regions). It takes every corner of the box bounded by the minimum and maximum of each feature, discards the corners that contradict a strong, significant correlation between two features, and projects the remaining corners with the PILOT matrix `A`. The convex hull of the projected corners is the estimated boundary.

Parts of the boundary far from any instance show where new test instances would extend the benchmark.

For a 2D projection the boundary is a closed polygon. For a 3D projection it is a triangulated convex surface; when the projected corners all lie on one plane (for example, a 3D projection of only two features), it is the flat polygon, triangulated.

## Examples

### Plot the boundary of a 2D instance space

```matlab
rootdir = 'test/data/example/';
if ~isfolder(rootdir), mkdir(rootdir); end
copyfile('test/data/metadata.csv', rootdir);
opts.perf = struct('MaxPerf', false, 'AbsPerf', true, 'epsilon', 0.20);
obj = InstanceSpace(rootdir, opts);
obj = obj.build('stages', {'prelim', 'sifted', 'pilot'});

out = CLOISTER(obj.model.data.X, obj.model.pilot.A, obj.opts.cloister);
Z = obj.model.pilot.Z;
plot(Z(:,1), Z(:,2), '.', out.Zecorr(:,1), out.Zecorr(:,2), 'r-')
legend('Instances', 'Estimated boundary')
```

### Plot the boundary of a 3D instance space

```matlab
obj.opts.pilot.dims = 3;
obj = obj.build('stages', {'pilot', 'cloister'});
c = obj.model.cloist;
trisurf(c.ZedgeFaces, c.Zedge(:,1), c.Zedge(:,2), c.Zedge(:,3), 'FaceAlpha', 0.1)
```

`obj.plot('boundary')` draws the same figure over the instances.

## Input Arguments

### `X` — Feature matrix

*numeric matrix*

Preprocessed, selected features, one row per instance. May contain `NaN`.

### `A` — Projection matrix

*`dims`-by-`nfeats` numeric matrix*

`model.pilot.A`.

### `opts` — Boundary options

*structure*

Normally `obj.opts.cloister`.

#### `opts.pval` — Significance level

*`0.05` (default) | scalar in (0, 1)*

Correlations with a larger p-value are ignored.

#### `opts.corrThreshold` — Correlation threshold

*`0.70` (default) | scalar in [0, 1]*

A corner is discarded when it contradicts the sign of a correlation stronger than this.

#### `opts.maxFeatures` — Feature limit

*`20` (default) | positive integer*

The number of corners grows as `2^nfeats`. With more features than this, CLOISTER warns (`ISA:CLOISTER:tooManyFeatures`) and returns the convex hull of the projected instances instead.

## Output Arguments

### `out` — Boundary

*structure*

#### `out.Zedge` — Boundary of all corners

2D: closed polygon, first vertex repeated last. 3D: vertices of the convex hull.

#### `out.Zecorr` — Boundary of the correlation-consistent corners

Same format as `Zedge`. Equals `Zedge` when no corner is discarded, when too many corners are discarded to form a hull, or when `opts.maxFeatures` is exceeded.

#### `out.ZedgeFaces`, `out.ZecorrFaces` — Hull triangulations

3D only: `nfaces`-by-3 indices into the rows of `Zedge` and `Zecorr`. Empty in 2D.

## Tips

- CLOISTER expects mean-centred features. `InstanceSpace` warns (`ISA:InstanceSpace:cloisterNotMeanCentred`) when `opts.auto.preproc` or `opts.norm.flag` is `false`.
- `scriptcsv` writes `Zedge` to `bounds.csv` and `Zecorr` to `bounds_prunned.csv`.
- If every projected corner lies on one line (a degenerate projection matrix or features), there is no region to bound and CLOISTER raises `ISA:CLOISTER:degenerateBoundary`.

## Version History

### v0.9.2 — 3D boundary

A 3D projection now gets a 3D convex hull, with `ZedgeFaces` and `ZecorrFaces`. Before, the hull used only the first two coordinates.

### v0.9.0 — Feature limit

`opts.maxFeatures` guards against an intractable number of corners.

## References

- Smith-Miles, K. & Muñoz, M.A. (2023). Instance Space Analysis for Algorithm Testing. *ACM Computing Surveys*, 55(12), Article 255. <https://doi.org/10.1145/3572895>

## See Also

`PILOT` | `InstanceSpace` | `scriptcsv`
