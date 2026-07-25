# Instance Space Analysis Toolkit — v0.9.0

This release is a complete refactor of the ISA Toolkit, covering everything since v0.3.3: a rewritten PILOT/PYTHIA/TRACE/SIFTED pipeline, a new class-based API, a full legacy-model migration path, and a repository-wide correctness and code-quality pass. It targets **MATLAB R2025a or later**.

---

## New functionality

**`InstanceSpace` class — the new primary API.** Wraps the pipeline in a single value class with:
- `build('stages', {...})` — run any subset of `{'prelim','sifted','pilot','cloister','pythia','trace'}`; prerequisites are checked and enforced automatically, in canonical order regardless of how they're listed.
- Options can be changed on the object between staged `build()` calls; re-running an earlier stage automatically invalidates every stage that transitively depends on it, so a save/explore/output step can never silently mix a freshly-recomputed stage with stale downstream results.
- `explore(testRootDir)`, `save()` / `InstanceSpace.load(rootdir)` (HDF5-compatible `-v7.3` `model.mat`), `getResults()`, `plot()`.
- `buildIS.m` and `exploreIS.m` remain as thin backward-compatibility wrappers around the class for existing callers (notably the MATILDA web platform).

**PILOT — 3D projection, alternative methods, and viewpoint optimisation.**
- `opts.pilot.dims = 3` for a 3D projection (analytic and numerical branches both generalised); `opts.pilot.dims = 2` remains the default.
- `opts.pilot.method = 'pls'` — Partial Least Squares as an alternative to the standard BFGS/analytic projection, useful when the feature matrix isn't full column rank.
- `opts.pilot.alpha` — adjustable weight on the performance-reconstruction term of PILOT's cost function.
- `PILOTviewpoint` — finds the optimal 2D camera viewpoint(s) of a 3D projection; `opts.pilot.viewGroups` requests one viewpoint per algorithm group instead of a single global one.

**PYTHIA — generic classifier registry.** Hardcoded SVM/LIBSVM replaced with a registry of MATLAB-native classifiers (`opts.pythia.classifier`: `'knn'`, `'svm'`, `'tree'`, `'nb'`, `'linear'`, `'ensemble'`), resolved via `ISAgetClassifierFcn`. Hyperparameter tuning via `opts.pythia.tuning`: scrambled Sobol sequence (default), MATLAB `bayesopt`, or pre-supplied `opts.pythia.params`.

**TRACE3 — unified footprint algorithm.** A single canonical implementation for both 2D and 3D instance spaces, using PYTHIA's predicted labels directly (mandatory PYTHIA→TRACE coupling). The pre-refactor DBSCAN + alpha-shape algorithm is retained as `TRACE_legacy.m`, selectable via `opts.trace.method = 'legacy'` (2D only).

**SIFTED — automated feature selection**, promoted from the former `SIFTED2` prototype: a correlation filter followed by correlation-clustering + a genetic algorithm that picks one representative feature per cluster, scored by PILOT's analytic branch.

**3D output.** `scriptpng` now writes a coloured z1/z2/z3 compass overlay on every 3D plot, and (for 3D projections, `opts.outputs.fig`) a `.fig` file alongside each footprint PNG for interactive rotation. `ISArecallView(fig, groupIdx)` snaps a reopened figure — even from a different MATLAB session — back to its optimised viewpoint.

**`ISAmigrateModel`** — migrates any pre-refactor `model.mat` to the current field layout in one call, file-based (`ISAmigrateModel(rootdir)`, with automatic backup) or in-memory. Covers the full legacy migration table: `opts` struct renames and merges, `opts.perf.MaxMin`→`MaxPerf`, the `model.data.bestPerformace` typo, `pythia.svm`/`.knn`→`.classifiers`, legacy hyperparameter field renames, missing `completedStages` inference, and — where the original training data is still available — **retraining** any classifiers still in the legacy LIBSVM struct format (which have no `predict()` method and can't simply be relabelled) using the current registry, and **recomputing** any pre-refactor DBSCAN+polyshape footprints fresh via TRACE3.

**`ISAvalidateOpts`** — validates every user-supplied `opts` field before defaults are filled in, raising a clear `ISA:ISAvalidateOpts:*` error on the first invalid value instead of letting a typo surface many stages later as a confusing crash.

**Density-based instance subsetting (`FILTER`)** now scales to datasets on the order of 20,000 instances via a KD-tree range search, instead of a dense O(n²) pairwise-distance matrix.

**`liveDemoIS.m`** — a new interactive, stage-by-stage walkthrough of the pipeline through the `InstanceSpace` class, replacing the previous `liveDemoIS.mlx` Live Script (which predated the class and several since-completed pipeline changes). Open it in MATLAB's Live Editor.

**AI-assisted analysis — Claude Code skill.** `.claude/skills/instance-space-analysis/` adds a general ISA methodology reference plus an operational reference for this repository's `options.json` schema and `InstanceSpace` class API, so an AI assistant using Claude Code against this repository can help run, interpret, and debug the pipeline. Scoped to this MATLAB toolkit only — it does not cover MATILDA's web interface or the Python `pyInstanceSpace` package, though its general methodology content may still help interpret results from those.

---

## Better engineering

- **Repository reorganised** into `core/`, `output/`, `utils/`, and `deprecated/` (plain subdirectories, not a `+ISA` package namespace — a deliberate choice to keep functions directly, standalone-callable and preserve MATILDA compatibility), with `startup.m` for interactive use and automatic path setup from any of the three entry points.
- **Test suite**: `test_integration.m` is now an exhaustive option-coverage regression suite (every classifier, tuning strategy, 2D/3D, PLS, viewpoint groups, staged `build()`/`explore()`/save-load round-trips through the class API, and the full `ISAmigrateModel` legacy-migration table) run against the reference dataset. `example.m` was split out as a separate, minimal getting-started demo.
- **Console output** standardised to a consistent `[STAGE] message` format throughout, with detailed per-trial/per-iteration output gated behind `opts.general.verbose`.
- **Named constants and centralised defaults**: every previously-hardcoded magic number is now a documented `opts` field with a sensible default (`ISAdefaults.m`), including `opts.cloister.corrThreshold`/`maxFeatures`, `opts.pythia.kFold`, `opts.prelim.iqrMultiplier`/`nanThreshold`, and more.
- **Performance**: `FILTER`'s pairwise-distance loop replaced with a KD-tree range search (see above); `SIFTED`'s correlation-selection loop and `PRELIM`'s tie-breaking loop vectorised; `SIFTED`'s GA fitness cache moved from an unsafe global `containers.Map` to a `persistent` variable, and its nested `parfor` (silently serialising under MATLAB's no-nested-parallelism rule) flattened to a plain loop.
- **Dead code removed**: `scriptdisc.m` and the unused `older_scripts/` directory deleted outright; several stale commented-out code blocks and a no-op self-assignment cleaned out of `core/`/`output/`.
- **Documentation**: README rewritten to match the pipeline's actual execution order (previously, PRELIM's and SIFTED's option sections were listed *after* PILOT/CLOISTER/PYTHIA/TRACE, even though they run first); every internal-only design-document citation (`spec §X.Y`, phase numbers) stripped from public docstrings; `Contents.m` and `CITATION.cff` added for MATLAB Central File Exchange / Zenodo / GitHub "Cite this repository" integration.
- **`.gitignore`** fixed so any file placed under `test/data/` is trackable (previously only two specific filenames were excepted from the blanket ignore rule).

**Breaking option-schema changes** (all pre-existing `model.mat`/`options.json` files can be brought forward automatically with `ISAmigrateModel`):
- `opts.pythia.uselibsvm` removed (no longer read); `opts.trace.usesim` removed (the PYTHIA/TRACE coupling is now unconditional); `opts.sifted.NTREES` removed (leftover from a since-replaced scoring step). `opts.pythia.useknn` is deprecated but still honoured as a fallback.
- `opts.cloister.cthres` → `opts.cloister.corrThreshold`; `opts.pythia.cvfolds` → `opts.pythia.kFold`; `opts.parallel.flag`/`.ncores` → `opts.general.parallel`/`.ncores`.

---

## Bug fixes

A representative sample of the correctness issues resolved across the refactor (see git history for the complete list):

- **Crashes fixed**: a `parfor` sliced-variable violation in PILOT that broke every build under parallel execution; PYTHIA's `sobolset` constructor call using an invalid `'Scramble'` name-value pair; a missing `end` that broke `scriptfcn.m`'s `axisLimits`; single-class/degenerate-label algorithms crashing PYTHIA instead of being handled explicitly; `ensurePool()` crashing when no parallel pool existed yet.
- **Silent correctness bugs fixed**: `explore()`/`evaluateTestSet` wasn't applying `opts.selvars.feats`/`opts.selvars.algos` before bound-clipping, and didn't validate that `metadata_test.csv`'s feature columns matched the trained model's in count *or order* — both could silently misalign or corrupt test-time normalisation. PYTHIA's eval mode truncated its outputs to the trained model's algorithm count instead of padding to the full reconciled count, so `explore()` crashed on any algorithm present in `metadata_test.csv` but absent from training (a case `evaluateTestSet`'s own reconciliation and TRACE's eval mode already supported); its outputs now pad correctly (`Yhat=false`, precision/recall/accuracy=`NaN`) for algorithms with no trained classifier. `InstanceSpace.load()` unconditionally marked every stage as completed regardless of what was actually in the saved model. Re-running an earlier `build()` stage left stale, no-longer-valid downstream results (e.g. PYTHIA/TRACE computed from an old PILOT projection) looking valid; the fix invalidates only the stages that actually depend on the one just re-run, not simply every stage listed after it.
- **Numerical/RNG fixes**: reproducibility issues in PYTHIA's Sobol/Bayes tuning under parallel execution; PRELIM's `min(Y(:))` computed before NaN exclusion; FILTER's uniformity metric now guards the degenerate case (fewer than two retained instances, or all retained instances coinciding) instead of silently producing `NaN`/`Inf`.
- **Output fixes**: PNG figures rendering with a black background in headless MATLAB sessions (R2025a+ figure `Theme` inheriting the desktop's dark-mode setting); spurious `InvertHardCopy` warnings; an interactive axes toolbar leaking into exported PNGs; missing gridlines and a fixed ±1 axis margin that crushed small-scale 3D projections (e.g. PILOT's PLS method) into a tiny fraction of the plot.
- **Migration-utility fixes** (found during this release's own review cycle): `ISAvalidateOpts` no longer silently bypasses validation for a field explicitly supplied as `[]` (distinct from "not supplied at all"); it now rejects `Inf`/`-Inf` consistently with the numeric checks downstream code already performs; `ISAmigrateModel`'s LIBSVM-classifier-name matching is now case-insensitive and accepts both `char` and `string` values.

---

## Change of licence and implications

**The repository's licence has changed from GPL v3 to the [PolyForm Noncommercial License 1.0.0](https://polyformproject.org/licenses/noncommercial/1.0.0/).** This is a breaking change for any existing fork or downstream user who was relying on GPL v3 terms.

What this means in practice:
- **Non-commercial use only, without prior permission.** PolyForm Noncommercial permits use, modification, and redistribution of the toolkit **for non-commercial research and educational purposes only**. GPL v3 placed no such restriction — commercial use of the previous release was permitted (subject to GPL's copyleft obligations). Any commercial use of this release or later requires prior written permission from the copyright holders.
- **No copyleft.** Unlike GPL v3, PolyForm Noncommercial does not require derivative works to be distributed under the same licence, or their source made available. You may build closed-source, non-commercial tools on top of this toolkit without a share-alike obligation.
- **Explicit patent and liability protections** are included in the new licence text that GPL v3 did not provide in the same form.
- **Every source file** now carries an `SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0` header, plus a one-sentence human-readable summary and the same copyright/reference block, replacing the previous bare GPL notices.
- The full licence text is at the repository root (`LICENSE`); `CITATION.cff` (new in this release) also records the licence for GitHub's "Cite this repository" integration and Zenodo.

**If you are already using a pre-v0.9.0 copy of this toolkit under GPL v3**, that copy remains under GPL v3 as licensed; the licence change applies to this release and all versions published from this point forward.
