# PR 64 Copilot follow-up

All four findings from the Balanced Copilot review were confirmed and addressed.

| Finding | Correction | Validation |
|---|---|---|
| Serial SIFTED dispatched worker cache resets | Resolve the pool once under the invocation's parallel setting, then pass only that pool to cache cleanup. Serial calls reset the client cache only. | Profile a serial call with an existing two-worker pool and verify no worker dispatch or pool deletion; rerun serial and parallel cache-isolation tests. |
| Legacy PLS boundaries remained uncentred | When recovering a missing PLS feature mean, translate saved CLOISTER vertices by `-Xmean*A'`. Preserve face indices and empty/absent boundaries. Existing means prevent repeat translation. | Test 2D and 3D models with populated, empty, and absent boundaries; check unchanged fitted projection, face indices, and migration idempotence. |
| CLOISTER omitted its fourth argument from help | Add four-argument syntax, row-vector shape, zero default, and fitted PLS centring semantics to MATLAB help and reference documentation. | Regenerate HTML and validate its source and links. |
| PYTHIA described all bad-class scores as probabilities | Define score semantics in the output reference and MATLAB help, including non-probability scores, unavailable values, and mixed CV columns. | Regenerate HTML and validate its source and links. |

## Verification

Before the code fixes, the new tests failed on the two populated-boundary cases and serial worker dispatch, reproducing both code findings. The empty/absent migration cases passed.

After the fixes, **12/12 targeted tests passed**, with zero failed or incomplete tests, in MATLAB R2026a Update 5. The run comprised all seven CopilotReviewTest cases, four existing SIFTED stage cases (including serial/parallel cache isolation), and the existing PLS projection-mean regression. Per-test results are recorded in `copilot-verification.json`.

Documentation generation is current (36 files); all 31 HTML pages passed link validation. `git diff --check` passed. The complete long-running integration suite was not repeated locally for this follow-up; GitHub CI runs on the pushed revision. No merge was performed.
