# TRACE

Compute or evaluate algorithm footprints in the instance space.

## Syntax

```
Training mode (7 args):
out = TRACE(Z, Ybin, Yhat, P, beta, algolabels, opts)

Evaluation mode (8 args):
out = TRACE(Z, Ybin, Yhat, P, beta, algolabels, opts, trainedTrace)
```

## Description

Builds footprints using TRACE3 (default) or the legacy DBSCAN algorithm.

Note: PYTHIA always runs before TRACE in this pipeline (mandatory coupling), so TRACE never trains its own KNN classifier -- that would be the fallback algorithm when no PYTHIA predictions are available, but it is superseded whenever PYTHIA has already produced predictions, which is always. If Yhat is empty/skipped, Zu = {yi=1} directly with a warning; there is no KNN-training fallback in this version.

## Input Arguments

| Argument | Description |
|---|---|
| `Z` | (`ninst` x `ndim`) projected instance coordinates (2D or 3D) |
| `Ybin` | (`ninst` x `nalgos`) true binary performance labels |
| `Yhat` | (`ninst` x `nalgos`) PYTHIA predicted labels; [] if unavailable |
| `P` | (`ninst` x 1) best-algorithm index per instance |
| `beta` | (`ninst` x 1) logical: instance is easy (has a good algorithm) |
| `algolabels` | (1 x `nalgos`) cell array of algorithm name strings |
| `opts` | struct with fields (see ISAdefaults for defaults): |
| &nbsp;&nbsp;`opts.method` | 'trace3' (default) or 'legacy' |
| &nbsp;&nbsp;`opts.PI` | minimum purity threshold (0.6) |
| &nbsp;&nbsp;`opts.minInstances` | minimum instances for a valid footprint (4) |
| &nbsp;&nbsp;`opts.minAreaFrac` | footprint must exceed this fraction of space (0.01) |
| &nbsp;&nbsp;`opts.contra` | contradiction removal; legacy only; default true when method='legacy', false otherwise |
| `trainedTrace` | (optional) trained trace struct from a prior TRACE call; when provided, activates evaluation mode. |

## Output Arguments

| Field | Description |
|---|---|
| `out` | struct with fields: |
| &nbsp;&nbsp;`out.space` | convex-hull space footprint (measure, density, ...) |
| &nbsp;&nbsp;`out.good{nalgos}` | good-performance footprints |
| &nbsp;&nbsp;`out.best{nalgos}` | best-algorithm footprints |
| &nbsp;&nbsp;`out.hard` | beta-hard footprint (~beta instances) |
| &nbsp;&nbsp;`out.summary` | (nalgos+1 x 11) cell array performance table |

## Legacy method (`opts.method = 'legacy'`)

Selectable via `opts.method='legacy'` above, implemented in `TRACE_legacy.m`. Uses DBSCAN clustering followed by boundary/polyshape construction per cluster, instead of TRACE3's default algorithm. Contradiction removal between best-algorithm footprints is applied when `useContra = true` (the default for legacy mode).

Inputs mirror TRACE's own training-mode signature above, with one confirmed difference: there is no `Yhat` argument -- legacy mode does not consult PYTHIA predictions at all. Outputs mirror the same `out` struct shape (`space`, `good`, `best`, `hard`), but each footprint struct's area is stored as `.area`, not `.measure` -- `TRACE.m` itself normalises this naming difference after calling `TRACE_legacy`, so callers of `TRACE.m` never see the discrepancy.

## References

- Smith-Miles, K. & Munoz, M.A. (2023). Instance Space Analysis for Algorithm Testing. ACM Computing Surveys, 55(12), Article 255. <https://doi.org/10.1145/3572895>
- Simpson, C., Munoz, M.A., Kandanaarachchi, S. & Campello, R.J.G.B. (2025). ISA3: A 3-dimensional expansion of Instance Space Analysis. Machine Learning, 114, 240. <https://doi.org/10.1007/s10994-025-06871-5>
- Munoz, M.A. & Smith-Miles, K. (2017). Performance analysis of continuous black-box optimization algorithms via footprints in instance space. Evolutionary Computation, 25(4), 529-554. <https://doi.org/10.1162/EVCO_a_00194>

## See Also

[PYTHIA](PYTHIA.html) | [InstanceSpace](InstanceSpace.html)
