# v0.9.2 release preparation

Prepared 17 October 2026 from `master` at `ed4bbc26b3381fd97e45d2b32ff93513f7d675c4`.
This document records release preparation; it does not mean the tag, GitHub
release, File Exchange update or Zenodo archive has been published.

## Scope and compatibility

The release collects merged work through PR #66: documentation and correctness
fixes, independent INIT, explicit persisted TRACE boundary tolerance, and the
bounded serial Octave workflow. See `RELEASE_NOTES.md` and the generated
`doc/html/WhatsNew.html` for changes and upgrade guidance. Rebuild existing
experiments to obtain corrected fitted results. Keep original models and
metadata until the new results have been checked.

MATLAB R2025a is the minimum supported release and CI target. Local validation
uses R2026a Update 5. The full Octave workflow is validated with 11.3.0,
Statistics 2.0.0 and Datatypes 1.5.0. See `octave-compatibility.md` for the exact
support boundary; portable archives do not imply cross-runtime model loading.

Closed issue #30 proposed optional canonical projection rotation, but that feature
is not implemented in current PILOT and is excluded from this release. Do not
use the closed issue state as evidence of implementation.

Issue #63 remains open: numerical PILOT restart ranking can amplify small
BLAS/thread-level differences. The release does not promise cross-environment
bitwise reproducibility. `opts.pilot.topoWeight` remains disabled. Classifier CV
summaries are conditional on fitted preprocessing, selection, projection and
tuning, rather than an outer-CV estimate of the whole pipeline.

## Verification

The release branch changes version metadata and documentation, including
`Contents.m`; it does not change runtime algorithms. Generated HTML is committed
so downloaded archives work with MATLAB Help without a Python installation.

Documentation build and internal link checks:

```sh
python -m pip install -r doc/requirements.txt
python doc/generate.py --check
python .github/scripts/check_doc_links.py doc/html
```

MATLAB, from the repository root:

```matlab
example
test_integration
```

The integration runner writes `coverage.xml` and `junit.xml` and fails on any
failed or incomplete case. Octave validation uses the commands in the
compatibility guide plus `runWorkflowEdges`, `runArchiveSmoke` and
`runPropertySmoke`. Run the 2D and 3D pipeline checks with PNG enabled and a Qt
display. On this machine the usable local display is `:101`; the inherited
remote display name rejected authorization. Do not interpret that initial
configuration failure as a supported-runtime test pass.

Local checks completed:

- Generated HTML matches all sources (37 generated files); link validation
  checks 32 pages with zero problems.
- Octave 11.3.0 with the pinned packages passes the portable preprocessing,
  analytic/numerical/PLS projection, KNN/CV/tuning, geometry, workflow edge,
  archive and property contracts. Both 2D and 3D full build/save/load/explore
  checks pass with CSV and PNG on the local Qt display. The initial remote
  display configuration failed its 3D graphics prerequisite; the Qt rerun
  passes without errors.
- Octave 6.4.0 passes the separate core-only analytic and preprocessing contracts;
  this is not full-pipeline support on that runtime.
- The source candidate archive passes ZIP integrity, version and required-file
  checks and excludes caches, environments and generated test/package outputs.
  Archive checksums and command logs are retained in ignored
  `.cache/release-0.9.2/` (runtime logs originate under `/tmp/`). Rebuild the
  candidate after changes and the final archive from the reviewed merge commit.

The full MATLAB result and CI status are recorded in
[release-preparation PR #67](https://github.com/andremun/InstanceSpace/pull/67).
Require a passing full MATLAB suite and green R2025a/pinned Octave CI checks on
the reviewed PR before publication.

## Publication checklist

1. Merge the release-preparation PR after required CI and Copilot review.
2. Confirm the final release date. If publication is later than 17 October 2026,
   update `Contents.m` and `CITATION.cff` consistently and regenerate/check the
   HTML. Their version is `0.9.2`; the release tag is `v0.9.2`.
3. Build a clean source archive from the reviewed commit, containing `startup.m`,
   `Contents.m`, `CITATION.cff`, `info.xml`, `LICENSE`, the pipeline code,
   generated `doc/html`, examples and tracked reference metadata. Exclude local
   environments, caches, generated test outputs and pinned package installs.
   Record its SHA-256 checksum and inspect its file list before uploading.
4. Create the `v0.9.2` tag on that exact commit and publish the GitHub release
   using only the v0.9.2 section of `RELEASE_NOTES.md`. State the rebuild guidance,
   MATLAB minimum, bounded Octave scope and known reproducibility limitation.
   Keep the existing release designation policy; preparation does not change
   the project's beta/stability status.
5. Verify the deployed reference site reports v0.9.2 and its Help-browser copy
   opens from the downloaded archive. Smoke-test the example in a clean folder.
6. Update MATLAB Central File Exchange and verify that Zenodo archives the new
   release. Retain the concept DOI; add a version-specific DOI only when issued.
   Verify metadata against the actual published artifacts.

Steps 1–6 are publication work, separate from preparing this branch.
