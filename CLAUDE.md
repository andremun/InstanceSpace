# ISA Toolkit — Claude Code Context

## Repository
- **Repo:** `andremun/InstanceSpace` (MATLAB, R2025a+)
- **Owner:** sole developer (Mario Andrés Muñoz) — no human reviewers; GitHub Copilot is the automated PR review gate. `master` has branch protection requiring changes go through a PR.
- **Primary API:** the `InstanceSpace` class (`InstanceSpace.m`) — `build()`/`explore()`/`save()`/`InstanceSpace.load()`/`getResults()`/`plot()`. `buildIS.m`/`exploreIS.m` are thin backward-compatibility wrappers around it (kept for callers like the MATILDA web platform); don't add new logic to them directly.
- **Completed refactor plan:** `docs/tech/isa_refactor_plan_v1.7.pdf` — historical source of truth for the v0.9.0 refactor (all ten phases). Now fully shipped; treat it as an archival reference, not a live backlog.
- **Current backlog:** GitHub Issues, milestone `v0.9.1` (issues #25–39 at time of writing) — this is the live source of truth for what's next, not the old plan document. See "Current state" below for the shape of it.
- **Reference papers:** `docs/papers/Smith-Miles2023.pdf` (general ISA methodology), `docs/papers/Simpson2025.pdf` (ISA3, the 3D extension — cited only in `PILOT.m`, `PILOTviewpoint.m`, `TRACE.m`, and the 3D-output files in `output/`, not repo-wide). PILOT's origin (Muñoz et al. 2018), legacy TRACE's origin (Muñoz & Smith-Miles 2017), and FILTER's origin (Alipour et al. 2023) are cited in each file's own docstring `Reference:` block but aren't archived as PDFs here.

## Test harness
Run from the repo root; `buildIS.m`/`exploreIS.m` print `EOF:SUCCESS` on success (errors abort before printing), and `test_integration.m` prints `EOF:SUCCESS`/`EOF:ERROR` at the end:
    example
    test_integration
`example.m` already runs the full `buildIS`+`exploreIS` pipeline internally on the bundled reference dataset — it's the one to start with. `test_integration.m` is the exhaustive option-coverage regression suite (every classifier, tuning strategy, 2D/3D, PLS, viewpoint groups, staged `build()`/`explore()`/save-load round-trips, the full `ISAmigrateModel` legacy-migration table); good as a reference for how an option is meant to be used, not the place to start reading.

## Current state (v0.9.0 shipped; v0.9.1 in progress)
- **All ten refactor-plan phases shipped as v0.9.0** (July 2026): rewritten PRELIM/SIFTED/PILOT/PYTHIA/TRACE pipeline, the `InstanceSpace` class API, `ISAmigrateModel`/`ISAvalidateOpts`, 3D projection/viewpoint support, and a licence change from GPL v3 to PolyForm Noncommercial 1.0.0. See `RELEASE_NOTES.md` for the full changelog.
- **v0.9.1 is next**: engineering-quality, architecture, and infrastructure follow-ups only — explicitly **no algorithmic changes, no new user-facing options**. Don't expand scope on one of these issues into a feature change without checking with the user first. Highlights:
  - Data ingestion gets a PYTHIA/TRACE-style train/eval structure via a new `INIT.m` (#38) — supersedes the narrower approaches originally proposed in #25 and #37; implement via #38, not those two separately.
  - CI (#34) and migrating `test_integration.m` to `matlab.unittest` (#39).
  - Two genuinely-still-open items from the old refactor plan's own "deferred to Phase 9" list, confirmed never actually closed when Phase 9 shipped: the `traceAlphaBoundary` multi-region alpha-shape bug (#31) and CLOISTER boundary display (#32).
  - `SECURITY.md`/`CONTRIBUTING.md`, LIBSVM MEX-file provenance, per-stage `build()`/`explore()` inspection callbacks, per-stage contract validation, PILOT projection-orientation canonicalisation, and a cross-repo reference-fixture export tool for the Python port (`pyInstanceSpace`) — see the milestone for the rest.
- **Repository layout**: `core/` (PRELIM, SIFTED, PILOT, PILOTviewpoint, CLOISTER, PYTHIA, TRACE, TRACE_legacy, FILTER), `output/` (scriptcsv, scriptpng, scriptweb, scriptfcn, ISArecallView), `utils/` (ISAdefaults, ISAvalidateOpts, ISAgetClassifierFcn, ISAmigrateModel, ISAsubsetData), `deprecated/` (PYTHIA2, PYTHIAtest, SIFTED2 — warn-and-forward shims kept for backward compatibility only; don't build new code on these).

## Workflow
- Work on feature branches; merge to `master` after Copilot review. `master`'s branch protection requires a PR — don't push directly to it.
- Fix root causes, not symptoms — no `real()` wrappers, no `isfinite()` guards.
- For v0.9.1-era work, cite the relevant GitHub issue number in commit messages (e.g. `#38`), not a refactor-plan section — the plan document is closed out; issues are the live backlog now.

## Key conventions
- Language: MATLAB; style follows standard MATLAB conventions.
- No magic numbers — use named `opts` fields with documented defaults (`ISAdefaults.m`), not hardcoded literals.
- Standalone pipeline functions (`PRELIM`, `SIFTED`, `PILOT`, `CLOISTER`, `PYTHIA`, `TRACE`, `FILTER`) must remain callable independently of the `InstanceSpace` class — deliberate, not an oversight, to preserve MATILDA web-platform compatibility. Don't fold them into class methods.
- **Train vs. eval mode is dispatched via `nargin`, not a fixed positional argument**: an optional trailing `trainedModel`/`trainedTrace`-style struct as the *last* argument selects eval mode over train mode (e.g. `PYTHIA(Z,Y,Ybin,Ybest,algolabels,opts)` trains; `PYTHIA(...,opts,trainedModel)` evaluates; same shape for `TRACE`). Only `PYTHIA` and `TRACE` currently have this. Confirmed by direct audit (see #38) that `SIFTED`/`PILOT`/`CLOISTER`/`FILTER` never recompute anything at explore time — they just reuse a stored trained artifact directly (e.g. `PILOT`'s trained projection is reused via a plain `X*A'` multiply) — so they don't need this structure and shouldn't get it added speculatively. `PRELIM` is the one confirmed gap (its explore-time equivalent is a separately-maintained, already-drifted reimplementation inside `InstanceSpace.evaluateTestSet`, not a shared function) — being fixed via `INIT.m` in #38.
- `opts.pilot.topoWeight` — still reserved/disabled, no effect in the current version; don't wire it up without confirming it's actually in scope for the task at hand.

## Known deferred items (tracked, not stale — check the linked issue before touching)
- `traceAlphaBoundary` (in `output/scriptcsv.m`) — silently exports partial boundaries for multi-region alpha shapes; its own comment admits it "works correctly for simple (single-region) connected alpha shapes" only. Was "deferred to Phase 9" in the old plan; Phase 9 shipped in v0.9.0 without this actually being fixed. Tracked as #31.
- CLOISTER's boundary is computed but never rendered in automated output (`scriptpng.m` has no reference to it at all; `InstanceSpace.plot()` has no `'boundary'` case). Also "deferred to Phase 9" in the old plan, also never actually closed. Tracked as #32.

Not still open, despite also appearing on the old plan's deferred list — closed by v0.9.0 but not previously recorded as such (#33, now documented in `RELEASE_NOTES.md`): the "ensemble hyperparameter API" question (verifying `fitcensemble`'s `OptimizeHyperparameters` support). PYTHIA's shipped design tunes every classifier, including `'ensemble'`, through one shared Sobol/Bayesian layer instead of `fitcensemble`'s own optimiser, sidestepping the question rather than answering it directly — don't re-open this as if it were still pending.

## `.gitignore` notes
- `docs/**/*.pdf` — tracked (negation rule), recursive under `docs/`.
- `test/data/*` — ignored except `metadata.csv`, `metadata_test.csv`, and `workspace_test.mat` (an explicit three-file allowlist). Do not broaden this to `!test/data/**` — that was tried and reverted because it makes every generated test-run output file (model.mat, CSVs, PNGs, per-case subdirectories from `test_integration.m`) trackable.
