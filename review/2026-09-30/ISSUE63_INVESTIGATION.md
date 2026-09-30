# Issue 63 investigation

**Conclusion: do not count PR #64 as fixing issue #63, and leave the issue open.** The reported large run-to-run differences were not reproduced on this host. The PR does not change the numerical solver or candidate-selection rules implicated by the issue. A passing local repetition is not evidence that the original CI problem has been removed.

## Source findings

Compared `master` (`08044f5`) with PR code at `f238f16`:

- PILOT's numerical objective function, `fminunc` loop and settings, starting-point generation, and restart-ranking block are unchanged. Ranking remains `[~,idx] = max(out.perf)`, where `perf` is a distance correlation rather than the reconstruction objective.
- PYTHIA still chooses the Sobol candidate with the minimum raw CV error. This PR changes failure handling and score metadata but introduces neither a near-tie tolerance nor a single-thread execution mode.
- The inspected `pyis_export_reference_data.m` requests `general.parallel=false`. With no pre-existing pool, master already uses zero workers. Respecting serial settings in this PR therefore does not explain or repair variation in a clean, pool-free export.
- On the bundled data with exporter defaults, SIFTED retains nine features and skips its clustering/GA because `K=10`. The persistent GA cache fix is not exercised by that case either.

## Controlled diagnostic

MATLAB R2026a Update 5, Linux GLNXA64. Prepare inputs once on master, using the reference metadata and the exporter's seed 42, serial execution, and default preprocessing/SIFTED options. This produces **212 instances, 9 selected features, and 10 algorithms**. Save the matrices in MAT format and pass those exact matrices to every PILOT call, with ten numerical restarts.

Run six fresh MATLAB sessions, with two PILOT calls in each session:

1. Master, default computation threads (12).
2. PR, default computation threads (12).
3. Master, `-singleCompThread` (1).
4. PR, `-singleCompThread` (1).
5. PR, another fresh default-thread session.
6. PR, another fresh single-thread session.

No parallel worker pool was present. This isolates PILOT; it is not a pair of full fixture exports. The existing downstream fixture bundle was not used as a numerical oracle because its manifest describes a different toolkit commit and macOS environment.

| Comparison | Result |
|---|---|
| Repeated calls in each of six sessions | Exact elementwise equality of `X0`, `alpha`, `A`, `Z`, `eoptim`, and `perf` |
| Master versus PR, matched thread setting | Exact equality of all six arrays |
| PR across fresh sessions, matched thread setting | Exact equality of all six arrays |
| Default versus single computation thread | `X0`, `alpha`, `A`, `Z`, and `eoptim` exactly equal; maximum absolute `perf` difference **6.5503158452884236e-15** |
| Selected restart | **10 in all 12 calls** |
| Gap between highest and second-highest restart score | About **0.00591624058314** in every call |

Thus thread count affects a few last bits of the correlation scores on this host, but not the solver coefficients or the selected restart. The measured score perturbation is many orders of magnitude below the winning gap. These observations do not establish the issue's proposed chain from BLAS scheduling to restart flips to SVM changes.

`issue63-verification.json` contains execution metadata, input hashes, winning gaps, and every comparison. `issue63_prepare_inputs.m` and `issue63_probe.m` reproduce the diagnostic. Use a separate MATLAB process per code/thread configuration; add this review directory to the MATLAB path before calling the helpers. For example:

```matlab
issue63_prepare_inputs('/path/to/master', '/tmp/isa63-inputs.mat')
% In a fresh session (with or without -singleCompThread):
issue63_probe('/path/to/toolkit', '/tmp/isa63-inputs.mat', '/tmp/isa63-pr-default')
```

## What remains unresolved

The original CI exports and runtime conditions must be reproduced before choosing a mitigation. Capture hashes of inputs entering each stage, every restart's `X0`/`alpha`/`perf`, the winning restart index, MATLAB update/BLAS version, computation-thread counts, and pool/worker settings. Compare both default-thread and single-thread full exports in that same environment. If PILOT inputs differ, trace preprocessing first; if only coefficients differ, investigate the solver; if the winner differs, measure the competing score gap before introducing a tolerance.

A single-thread run cannot be called a fix solely because it passes on a host where the default-thread run also passes. Likewise, `ntries=1` removes restart selection but does not guarantee numerical bit-reproducibility; fixed PYTHIA parameters remove tuning selection, not all numerical variation. No tolerance, solver setting, or performance definition was changed for this investigation.
