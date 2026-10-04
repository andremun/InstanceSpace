# GNU Octave serial release

Updated 4 October 2026. The supplied licensing/Octave plan was used as an
architecture proposal, not as an instruction to change licensing or distribute
third-party software. The implementation keeps the MATLAB backends and places
Octave-specific algorithms and runtime operations in `utils/+isacompat`.

## Supported workflow

The release scope is a serial **2D/3D build → save/load → explore** workflow:

| Stage | Octave implementation and boundary |
| --- | --- |
| INIT | Datatypes CSV tables, instance/feature/algorithm labels and categorical sources; existing train/evaluation schema validation. |
| PRELIM | Missing-data rules, bounds, performance labels, fitted Box–Cox maximum likelihood and application of frozen normalization parameters. Constant positive Box–Cox input uses lambda 1 because the likelihood is unidentifiable. |
| SIFTED | Correlation screening, correlation-distance k-means, CV/KNN fitness, exact enumeration of small integer spaces or a serial integer GA. Clustering uses batch updates; MATLAB retains its online phase. The selected search backend and fitness are recorded in `sifted.search`, and partitions in `sifted.cvpartition`. |
| PILOT | Analytic and numerical standard projection, SIMPLS projection, 2D/3D and optimized 3D camera views. Octave uses serial `fminunc` with mapped options. |
| CLOISTER | Pearson significance and full/correlation-constrained hulls, including 3D facets and degenerate cases. |
| PYTHIA | **KNN** is the validated classifier: weighted/unweighted fitting, constant labels, CV, `none` and `sobol` tuning, held-out evaluation and missing-outcome scoring. Other classifier names fail explicitly. |
| TRACE | TRACE3 regularized Delaunay alpha complexes in 2D/3D, containment, measures, boundaries, disconnected regions, holes and region-size thresholds. Legacy TRACE is unsupported. |
| Persistence | Versioned Octave archive with explicit reconstruction; classifier archives remain runtime-specific. MATLAB retains its existing MAT format. |
| Output | CSV and PNG; Qt is required for 3D PNG. Octave defaults `outputs.fig=false`. |

Bayesian tuning, parallel execution, `.fig`, Live Editor, the web palette export,
legacy polygon operations and unrestricted cross-runtime classifier loading are
outside this release. Unsupported pipeline options fail with an actionable
compatibility error. There is no automatic dependency installation or loading
in toolkit startup/build.

## Reproduce the tested environment

The full pipeline was exercised with **Octave 11.3.0**, **Statistics 2.0.0** and
**Datatypes 1.5.0**. The minimum package runtime is Octave 11.1; this does not
claim that every 11.x version has been tested. The tested Flatpak application
commit is `a37e6c8b9898ad580b4666407cb324fe8cd6e18dd1860fa811028b34575c577f`
with the KDE 6.10 runtime. The isolated test helper pins package URLs, versions
and SHA-256 checksums. Package archives/builds stay in ignored test directories.

From the repository root:

```sh
# Explicit opt-in installation of the pinned test dependencies.
flatpak run org.octave.Octave --no-gui --quiet --eval \
  "addpath('tests/portable'); setupOctaveValidation(true);"
# Subsequent sessions load the isolated environment, without reinstalling.
flatpak run org.octave.Octave --no-gui --quiet --eval \
  "addpath('tests/portable'); setupOctaveValidation(); runPortableSmoke(); runLearningSmoke(); runOptimizationSmoke(); runGeometrySmoke();"
```

For a native installation, replace `flatpak run org.octave.Octave` with its
`octave` executable. `/usr/bin/octave` may still be the older distribution
runtime. The test helper selects `test/data/octave-validation/octave_packages`
for that session; an ordinary session does not automatically load it.

A complete application example after loading the packages:

```matlab
opts.general.parallel = false;
opts.pythia.classifier = 'knn';
opts.pythia.tuning = 'sobol';
opts.pythia.nTuningIter = 20;
opts.pilot.dims = 3;                    % or 2
opts.outputs.fig = false;
obj = InstanceSpace('path/to/training', opts);
obj = obj.build();                      % metadata.csv; writes model.mat
obj = InstanceSpace.load('path/to/training');
obj = obj.explore('path/to/test');       % metadata_test.csv
```

Use a desktop display or a virtual display for Qt rendering. For headless CI:

```sh
xvfb-run -a flatpak run org.octave.Octave --no-gui --quiet --eval \
  "addpath('tests/portable'); setupOctaveValidation(); set(0,'defaultfigurevisible','off'); runPipelineSmoke(2,true); runPipelineSmoke(3,true);"
```

`scriptpng` selects Qt when available. It rejects 3D export with gnuplot because
that backend produced invalid or unusable footprint/legend plots in testing.
Octave uses the viridis color map; MATLAB retains parula. PNG pixel identity is
not a cross-runtime contract.

## Numerical and geometry contracts

Equal random seeds do **not** guarantee equal folds, selected subsets, optimizer
trajectories, tie choices or classifiers between MATLAB and Octave. Both preserve
caller RNG state where the existing stage API promises it. Saved CV partitions
retain their masks; PYTHIA records the actual Sobol hyperparameter candidates in
`pythia.tuningCandidates`.

Octave's two-dimensional Sobol design uses Gray-code direction numbers, skips
index zero and applies a random digital shift. MATLAB continues to use its
Matousek–Affine–Owen scramble. They are different randomizations of a Sobol
sequence, not an assertion of identical candidates for equal seeds.

The Octave footprint consists of full-dimensional Delaunay simplices with
circumradius at most alpha. The default alpha is the smallest radius at which
every input point belongs to a retained full-dimensional simplex. Regions join
across full facets, holes are retained, region thresholds remove components by
area/volume, and containment includes boundaries. Singular edges/faces do not
contribute. Duplicate points are removed; rank-deficient supports produce empty
footprints. Cospherical inputs use Qhull's triangulation with the point-at-infinity
option. Boundary orientation/triangulation and roundoff choices can differ from
MATLAB `alphaShape`; exact arbitrary-data geometry parity is not claimed.

Summary rounding, NaN-aware extrema, midpoint quartiles and Pearson p-values
have explicit adapters. Box–Cox uses a centered log-likelihood to avoid scale
bias and `expm1` for a stable transformation near lambda zero. SIMPLS uses centered
scores and loadings with `Z=(X-Xmean)*A'`.

## Archive contract

Octave writes MAT v7 binary data with `archiveVersion=1` and an encoded `payload`.
Arrays, cell arrays and structs remain data; categories, CV masks and alpha
complex inputs have explicit type markers and reconstruction paths. KNN models
use the Statistics package's `savemodel`/`loadmodel` API, embedded as byte arrays:
these are **Octave/Statistics-specific classifier archives**, not portable MATLAB
classifier objects. Keep the pinned runtime/packages to reload them. Reconstruction
does not retrain classifiers. Unknown types/schema versions fail explicitly.
Writes use a temporary file followed by replacement after serialization succeeds.
`InstanceSpace.load` and file-based `ISAmigrateModel` understand the archive.

MATLAB keeps flattened top-level model fields and `-v7.3` persistence. Exported
CSV numerical results are the supported interchange surface; do not load an
Octave classifier archive as a MATLAB model.

## Validation and limits

- `runPortableSmoke`: independent SVD reconstruction, weights, R-squared,
  summaries, PRELIM bounds/normalization/labels and CLOISTER 2D/3D.
- `runLearningSmoke`: weighted and unweighted KNN, class order, disjoint/exhaustive
  CV masks, probabilities, constant labels, missing test outcomes and Sobol tuning.
- `runOptimizationSmoke`: centered PLS scores, numerical projection consistency,
  caller RNG, bounded integer search and Sobol coverage/repeatability.
- `runGeometrySmoke`: square/cube/simplex measures, containment, duplicate and
  rank-deficient points, disconnected regions, thresholds and an annular hole.
- `runPipelineSmoke`: quoted Unicode instance IDs, sources, normalization, feature
  selection, full 2D/3D build, save/load/explore, matching predictions/summaries,
  evaluation-only algorithms, CSV and PNG signatures. Output lives in ignored test directories.
- `runWorkflowEdges`: fractional/density subsets and preservation of the previous
  archive when serialization fails.

The shared deterministic/learning/optimization tests also passed on MATLAB
R2026a Update 5. Both existing `ClassApiTest` cases passed, including staged
build/invalidation, save/load/explore and callbacks. They were completed in
separate runs after a combined run hit its graphics-output time limit.
The existing MATLAB CI remains the check for supported R2025a. A separate minimal
Octave 6.4 job exercises core-only deterministic contracts; it is not a supported
full-pipeline runtime. The modern Octave CI job pins the Flatpak application and
package versions and runs the serial suites and complete PNG workflows.

This is a bounded serial release, not a claim of every MATLAB option or arbitrary
cross-runtime numerical parity. No historical reference fixtures were replaced.
Statistics/Datatypes are separately obtained dependencies; their code or binaries
are not bundled here and no toolkit licence has changed.
