# PILOT

Project instance features onto a 2D or 3D instance space

<!-- opts: pilot -->

## Syntax

```
out = PILOT(X,Y,featlabels,opts)
```

## Description

`out = PILOT(X,Y,featlabels,opts)` finds a linear projection `Z = X*A'` of the feature matrix `X` into `opts.dims` dimensions, chosen so that both the features and the algorithm performance `Y` can be reconstructed from `Z` as well as possible (Prediction-based Linear Dimensionality Reduction). Instances with similar features and similar performance land close together, so trends in performance become visible as directions in the space.

PILOT has three solvers, selected by `opts.method` and `opts.analytic`:

- **Numerical** (default): BFGS minimisation of the reconstruction error of `X` and `Y` from `Z` (see Algorithms), restarted `opts.ntries` times. The trial whose projection best preserves the pairwise distances of the original data is kept.
- **Analytic**: the closed-form eigenvector solution. Fast and deterministic. If `X` is rank deficient, PILOT falls back to the numerical solver.
- **PLS**: Partial Least Squares regression of `Y` on `X` with `plsregress`. `opts.alpha` does not apply.

## Examples

### Project the reference data in 2D

Prepare the features with the earlier pipeline stages, then call PILOT directly.

```matlab
rootdir = 'test/data/example/';
if ~isfolder(rootdir), mkdir(rootdir); end
copyfile('test/data/metadata.csv', rootdir);
opts.perf = struct('MaxPerf', false, 'AbsPerf', true, 'epsilon', 0.20);
obj = InstanceSpace(rootdir, opts);
obj = obj.build('stages', {'prelim', 'sifted'});

d = obj.model.data;
out = PILOT(d.X, d.Y, d.featlabels, obj.opts.pilot);
scatter(out.Z(:,1), out.Z(:,2), 12, 'filled')
xlabel('z_1'), ylabel('z_2')
```

### Compare the solvers

The analytic solver is deterministic and much faster than BFGS. Compare how well each reconstructs the data.

```matlab
pilotOpts = obj.opts.pilot;
pilotOpts.analytic = true;
outA = PILOT(d.X, d.Y, d.featlabels, pilotOpts);

pilotOpts.analytic = false;
pilotOpts.method = 'pls';
outP = PILOT(d.X, d.Y, d.featlabels, pilotOpts);

fprintf('Analytic error: %.3f   PLS error: %.3f\n', outA.error, outP.error)
```

### Build a 3D instance space

Set `opts.pilot.dims` to 3 before building. `InstanceSpace` then also calls `PILOTviewpoint` to find the best 2D viewing angles.

```matlab
obj.opts.pilot.dims = 3;
obj = obj.build('stages', {'pilot'});
Z = obj.model.pilot.Z;
scatter3(Z(:,1), Z(:,2), Z(:,3), 12, 'filled')
```

## Input Arguments

### `X` — Feature matrix

*numeric matrix*

Instance features, one row per instance and one column per feature, normally the output of `PRELIM` and `SIFTED`.

### `Y` — Performance matrix

*numeric matrix*

Algorithm performance, one row per instance and one column per algorithm. Must have the same number of rows as `X`.

### `featlabels` — Feature names

*cell array of character vectors*

One name per column of `X`. Used for the row labels of `out.summary`.

### `opts` — Projection options

*structure*

Normally `obj.opts.pilot`. Missing fields take the defaults set by `ISAdefaults`.

#### `opts.dims` — Number of dimensions

*`2` (default) | `3`*

Dimension of the instance space.

#### `opts.method` — Projection method

*`'standard'` (default) | `'pls'`*

`'standard'` uses the analytic or numerical PILOT solution. `'pls'` uses Partial Least Squares.

#### `opts.analytic` — Use the closed-form solution

*`false` (default) | `true`*

Only used when `opts.method` is `'standard'`. Falls back to the numerical solver when `X` is rank deficient.

#### `opts.ntries` — Number of BFGS restarts

*`10` (default) | positive integer*

More restarts reduce the chance of a poor local minimum and take proportionally longer.

#### `opts.alpha` — Weight of the performance term

*`1.0` (default) | positive scalar*

Weight of `||Y - C*Z||^2` relative to the feature term. Larger values favour a projection that explains performance over one that explains features. Not used by `'pls'`.

#### `opts.seed` — Random seed

*`opts.general.seed` (default) | integer*

Seeds the BFGS starting points, so repeated runs give the same projection.

#### `opts.verbose` — Report progress

*`opts.general.verbose` (default) | logical*

#### `opts.X0` — Starting points

*numeric matrix*

Optional user-supplied BFGS starting points, one column per trial.

#### `opts.precalcAlpha` — Precomputed solution

*numeric vector*

Optional full BFGS solution vector from an earlier run. Skips the optimisation.

## Output Arguments

### `out` — Projection

*structure*

#### `out.A` — Projection matrix

`dims`-by-`nfeats` matrix. Project new, identically preprocessed instances with `Z = Xnew*out.A'`.

#### `out.Z` — Instance coordinates

`ninst`-by-`dims` coordinates of each instance in the instance space.

#### `out.B`, `out.C` — Reconstruction matrices

Reconstruct the features and the performance from `Z`.

#### `out.error` — Reconstruction error

Sum of squared errors of the reconstruction of `[X Y]`.

#### `out.R2` — Coefficient of determination

`R^2` of the reconstruction of each column of `[X Y]`.

#### `out.summary` — Projection table

Cell array showing `A` with the feature labels. Also written to `projection_matrix.csv` by `scriptcsv`.

#### `out.alpha`, `out.X0`, `out.eoptim`, `out.perf` — Trial results

Numerical solver only: the solution, start point, cost, and distance-preservation score of every trial.

## Algorithms

For the standard method, PILOT solves

```
min  ||X - Z*B'||^2 + alpha*||Y - Z*C||^2    subject to  Z = X*A'
```

where `X` and `Y` are the preprocessed feature and performance matrices. The analytic solution takes the leading `dims` eigenvectors of `[X Y]'*[X Y]`, with the `Y` block weighted by `sqrt(alpha)` (Muñoz et al., 2018). The numerical solution runs BFGS (`fminunc`) from `opts.ntries` random starts and keeps the trial with the highest correlation between the pairwise distances of the instances before and after projection, so the chosen space also preserves the topology of the data.

## Version History

### v0.9.0 — 3D projections and PLS

`opts.dims = 3` produces a 3D instance space (ISA3). `opts.method = 'pls'` adds Partial Least Squares as an alternative projection.

### v0.9.1 — Seed control

The BFGS starting points use `opts.seed`, so a changed `opts.general.seed` changes the projection and the same seed reproduces it.

## References

- Muñoz, M.A., Villanova, L., Baatar, D. & Smith-Miles, K. (2018). Instance spaces for machine learning classification. *Machine Learning*, 107(1), 109–147. <https://doi.org/10.1007/s10994-017-5629-5>
- Smith-Miles, K. & Muñoz, M.A. (2023). Instance Space Analysis for Algorithm Testing. *ACM Computing Surveys*, 55(12), Article 255. <https://doi.org/10.1145/3572895>
- Simpson, C., Muñoz, M.A., Kandanaarachchi, S. & Campello, R.J.G.B. (2025). ISA3: A 3-dimensional expansion of Instance Space Analysis. *Machine Learning*, 114, 240. <https://doi.org/10.1007/s10994-025-06871-5>

## See Also

`SIFTED` | `PILOTviewpoint` | `CLOISTER` | `InstanceSpace`
