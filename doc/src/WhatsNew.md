# WhatsNew

This page provides a condensed summary of recent changes and updates to the software. For full details, including bug fixes and the rationale behind the license change, please refer to `RELEASE_NOTES.md` in the repository root.

## v0.9.1 (Current Release)

- **Target Compatibility:** This release targets MATLAB R2025a or later.
- **Scope:** This is an engineering/infrastructure follow-up to v0.9.0, with no changes to pipeline algorithm behavior.
- **Headline Bullets:**
  - Introduced the `onStage` name-value argument for per-stage inspection callbacks in `build()`/`explore()`.
  - Implemented CLOISTER boundary rendering for 2D projections (`scriptpng.m`, `InstanceSpace.plot('boundary')`).
  - Engineering improvements include dual-mode data ingestion (`INIT.m`), per-stage I/O contract validation, GitHub Actions CI, and migration to `matlab.unittest`.
  - Key bug fixes include explore-time tie-breaking, multi-region footprint export, and seed propagation to PILOT/SIFTED.

## v0.9.0 (Previous Release)

- **Scope:** This release represents a complete toolkit refactor.
- **Headline Bullets:**
  - Introduced the `InstanceSpace` class as the primary API, replacing direct script calls with `buildIS.m`/`exploreIS.m` wrappers.
  - Major updates to PILOT (3D projection/PLS), PYTHIA (generic classifier registry), and TRACE3 (unified footprint algorithm).
  - New features include SIFTED, `ISAmigrateModel`, `ISAvalidateOpts`, and `liveDemoIS.m`.
  - Breaking changes include a license change to PolyForm Noncommercial 1.0.0 and option-schema updates (migratable via `ISAmigrateModel`).

See `RELEASE_NOTES.md` in the repository root for the complete list of changes, bug fixes, and the full licence-change rationale.
