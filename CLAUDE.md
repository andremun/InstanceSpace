# ISA Toolkit — Claude Code Context

## Repository
- **Repo:** `andremun/InstanceSpace` (MATLAB, R2025a+)
- **Owner:** sole developer (Mario Andrés Muñoz) — no human reviewers; GitHub Copilot is the automated PR review gate. `master` has branch protection requiring changes go through a PR.
- **Primary API:** the `InstanceSpace` class (`InstanceSpace.m`) — `build()`/`explore()`/`save()`/`InstanceSpace.load()`/`getResults()`/`plot()`. `buildIS.m`/`exploreIS.m` are thin backward-compatibility wrappers around it (kept for callers like the MATILDA web platform); don't add new logic to them directly.
- **Completed refactor plan:** `docs/tech/isa_refactor_plan_v1.7.pdf` — historical source of truth for the v0.9.0 refactor (all ten phases). Now fully shipped; treat it as an archival reference, not a live backlog.
- **Current backlog:** GitHub Issues are the live source of truth. As of 17 October 2026, #63 (PILOT numerical restart sensitivity to floating-point jitter) is the only open MATLAB issue. Do not treat old release plans as pending work.
- **Reference papers:** `docs/papers/Smith-Miles2023.pdf` (general ISA methodology), `docs/papers/Simpson2025.pdf` (ISA3, the 3D extension — cited only in `PILOT.m`, `PILOTviewpoint.m`, `TRACE.m`, and the 3D-output files in `output/`, not repo-wide). PILOT's origin (Muñoz et al. 2018), legacy TRACE's origin (Muñoz & Smith-Miles 2017), and FILTER's origin (Alipour et al. 2023) are cited in each file's own docstring `Reference:` block but aren't archived as PDFs here.

## Test harness
Run from the repo root; `buildIS.m`/`exploreIS.m` print `EOF:SUCCESS` on success (errors abort before printing), and `test_integration.m` prints `EOF:SUCCESS`/`EOF:ERROR` at the end:
    example
    test_integration
`example.m` already runs the full `buildIS`+`exploreIS` pipeline internally on the bundled reference dataset — it's the one to start with. `test_integration.m` is a thin `matlab.unittest`-based runner (#39) over `tests/*.m` (`PipelineOptionsTest` — every classifier, tuning strategy, 2D/3D, PLS, viewpoint groups, parameterised; `ClassApiTest` — staged `build()`/`explore()`/save-load round-trips; `MigrationTest` — the full `ISAmigrateModel` legacy-migration table; `RegressionTest` — targeted fixed-bug regressions), with a `CodeCoveragePlugin` producing `coverage.xml`. New features and bug fixes get a test *method* under `tests/`, not a new block in `test_integration.m` itself. Good as a reference for how an option is meant to be used, not the place to start reading. Runs automatically on every push/PR via CI (#34, `.github/workflows/tests.yml`); pinned to `release: R2025a`, not `latest` (deliberate — this is the documented minimum supported MATLAB release).

## Current state (v0.9.2 release preparation)
- v0.9.0 and v0.9.1 are published. The v0.9.2 changes are merged to `master` through #60, #62, #64, #65 and #66. Release metadata and documentation are being prepared separately from publishing a tag/release; use `docs/release-0.9.2.md` for the release checklist.
- Implemented and closed: #50 (3D CLOISTER), #52 (hole boundary cycles), #53 (reference site), #57 (JUnit/Codecov reporting), #58 (unobserved test outcomes), #59 (Oracle probability), #61 (SVM calibration seeding). These are not pending merges. #30 is also closed on GitHub, but its proposed optional canonical rotation is not implemented in current PILOT; do not infer implementation from the issue state. `RELEASE_NOTES.md` is the complete user-facing change record.
- Current contracts include independent INIT snapshots, preprocessing recomputation after portfolio pruning, SIFTED cache isolation and replay inputs, fitted PLS centering, training-only selection policy, honest CV/test summaries, stage provenance, persisted TRACE distance tolerance and serial Octave support. Review `doc/src/`, `docs/octave-compatibility.md` and `review/2026-09-30/FINAL_REPORT.md` before changing them.
- Repository layout: `core/` stages, `output/` exporters, `utils/` utilities and `+isacompat` runtime adapters, `deprecated/` warn-and-forward shims, `tests/` MATLAB tests and `tests/portable/` shared contracts. `.claude/skills/instance-space-analysis/` contains methodology and repository usage guidance. The older referenced `matlab-ci-testing` skill is not present in this checkout.
- CI targets MATLAB R2025a; local release validation uses R2026a Update 5. The full serial Octave workflow is validated with 11.3.0, Statistics 2.0.0 and Datatypes 1.5.0. Core-only checks on older Octave do not imply full workflow support.

## Workflow
- Work on feature branches; merge to `master` after Copilot review. `master`'s branch protection requires a PR — don't push directly to it.
- Fix root causes, not symptoms — no `real()` wrappers, no `isfinite()` guards.
- Cite the relevant GitHub issue number in commit messages (e.g. `#38`) for issue-tracked work, not a refactor-plan section — the plan document is closed out; issues are the live backlog.
- A docs-only change (no `.m` files) doesn't need the same review scrutiny as a code change — #54/#55 both merged directly once CI was green, no extended review cycle.

## Key conventions
- Language: MATLAB; style follows standard MATLAB conventions.
- No magic numbers — use named `opts` fields with documented defaults (`ISAdefaults.m`), not hardcoded literals.
- Standalone pipeline functions (`INIT`, `PRELIM`, `SIFTED`, `PILOT`, `CLOISTER`, `PYTHIA`, `TRACE`, `FILTER`) must remain callable independently of the `InstanceSpace` class — deliberate, not an oversight, to preserve MATILDA web-platform compatibility. Don't fold them into class methods.
- **Train vs. eval mode is dispatched via `nargin`, not a fixed positional argument**: an optional trailing `trainedModel`/`trainedTrace`/`trainedPrelim`-style struct as the *last* argument selects eval mode over train mode (e.g. `PYTHIA(Z,Y,Ybin,Ybest,algolabels,opts)` trains; `PYTHIA(...,opts,trainedModel)` evaluates; same shape for `TRACE` and, since #38, `PRELIM`/`INIT`). Confirmed by direct audit (see #38) that `SIFTED`/`PILOT`/`CLOISTER`/`FILTER` never recompute anything at explore time — they just reuse a stored trained artifact directly (standard PILOT reuses `X*A'`; PLS applies its stored feature mean before that multiply) — so they don't need this structure and shouldn't get it added speculatively.
- **Data ingestion (#38)**: `INIT.m` (`core/`) does the CSV-read + `opts.selvars.feats`/`.algos` filtering that used to be duplicated between `InstanceSpace.runPrelim` and `InstanceSpace.evaluateTestSet` — `INIT(rootdir, opts)` trains (reads `metadata.csv`), `INIT(rootdir, opts, trainedModel)` evaluates (reads `metadata_test.csv`, validates the feature set/order against `trainedModel`, reconciles algorithm columns). `PRELIM.m` itself now also has the `nargin` dual mode: its Ybin/Ybest/P/beta computation (including tie-breaking) is unconditional shared code, only the bound-clipping/Box-Cox+Z-score normalisation branches on train (fit) vs eval (apply trainedPrelim's already-fit parameters).
- `opts.pilot.topoWeight` — still reserved/disabled, no effect in the current version; don't wire it up without confirming it's actually in scope for the task at hand.

## Documentation conventions
- **`README.md` is user-facing only.** No GitHub issue numbers, no "confirmed via CI", no "an earlier version did X", no internal spec-section references (`spec §N`) — state the current fact, not how it got that way. That kind of context belongs in `RELEASE_NOTES.md`, `CLAUDE.md`, or commit messages/issues instead. (This was found drifting in during v0.9.1 and cleaned up in #55 — don't let it drift back in.)
- **Inline code in Markdown docs uses single backticks** (`` `like this` ``), not triple backticks. `README.md` had drifted to ~208 triple-backtick inline spans before #55's cleanup; `CONTRIBUTING.md`/`SECURITY.md`/`CLAUDE.md`/`RELEASE_NOTES.md` were already consistent. Triple backticks are still correct for actual fenced code blocks.
- **Every `.m` file's docstring** starts with the same copyright/SPDX-license header block (copy an existing file's, e.g. `core/PILOT.m`) followed by a `Reference:` section citing only the paper(s) that specific file implements. `core/` pipeline-stage files (`PRELIM`/`SIFTED`/`PILOT`/`PILOTviewpoint`/`CLOISTER`/`PYTHIA`/`TRACE`/`TRACE_legacy`/`FILTER`/`INIT`) and `InstanceSpace.m` additionally use a structured `Inputs`/`Outputs` section (see any of them for the exact format); `utils/`/`output/`/`deprecated/`/`buildIS.m`/`exploreIS.m` intentionally use lighter narrative-prose docstrings instead — that split (algorithm implementations vs. supporting utilities) is deliberate, not an inconsistency to "fix."

## Known limitations
- **#63 remains open:** numerical PILOT restart selection can amplify BLAS/thread-level floating-point jitter. Fixed-seed fixture repeats in a pinned environment do not establish cross-environment reproducibility or resolve this issue.
- `opts.pilot.topoWeight` remains disabled. Do not silently activate it as documentation cleanup.
- Octave supports only the bounded serial workflow documented in `docs/octave-compatibility.md`; its model archives are runtime-specific. The TRACE boundary policy agrees with Python, but arbitrary numerical/geometry identity is not promised.

## Documentation site (#53)
- Sources: `doc/src/*.md`; generator: `doc/generate.py` (python-markdown, pinned in `doc/requirements.txt`); output: `doc/html/` (committed — the MATLAB Help browser reads it via `info.xml` from a plain clone).
- The `TOC` list in `doc/generate.py` is the single source for the sidebar, breadcrumbs, `helptoc.xml`, the generated *All Functions* page and search. A page not in `TOC` fails the build.
- Page conventions (Syntax/Description/Examples/Input Arguments/Output Arguments/Version History/See Also, ``### `name` — Short description`` argument headings, `<!-- opts: pilot -->` for bare `opts.x` links) are documented in the generator's docstring.
- `.github/workflows/docs-pages.yml`: `generate.py --check` + `.github/scripts/check_doc_links.py` on PRs; deploys `doc/html` to GitHub Pages from `master`. Deployment uses the repository's GitHub Pages settings (Source: GitHub Actions).
- When a function's behaviour or signature changes, update its `doc/src` page (and a Version History entry), then regenerate.

**LIBSVM MEX-files removed (#29).** `svmpredict.mexw64`/`svmtrain.mexw64` were precompiled binaries with no corresponding source in the tree — removed entirely rather than kept-with-provenance, since LIBSVM is already fully deprecated for new runs. `PYTHIA`'s eval mode now raises a clear `ISA:PYTHIA:noLibsvm` error (naming the algorithm and pointing to the official LIBSVM project) if it ever needs to dispatch to a legacy LIBSVM-format classifier struct and `svmpredict` isn't on the path, instead of MATLAB's generic undefined-function error. Don't re-add these binaries; if LIBSVM support is ever needed again, get it from upstream rather than re-committing a precompiled MEX-file with no visible build source.

## Cross-repo note
`andremun/pyInstanceSpace` is the Python implementation. Keep Python-specific changes and fixture tooling there. Python issue #367 tracks improving its documentation website toward this reference site's navigation and presentation. Shared boundary policy does not imply full numerical parity. Check current connector permissions rather than relying on historical access failures.

## `.gitignore` notes
- `docs/**/*.pdf` — tracked (negation rule), recursive under `docs/`.
- `test/data/*` — ignored except `metadata.csv`, `metadata_test.csv`, and `workspace_test.mat` (an explicit three-file allowlist). Do not broaden this to `!test/data/**` — that was tried and reverted because it makes every generated test-run output file (model.mat, CSVs, PNGs, per-case subdirectories from `test_integration.m`) trackable.
