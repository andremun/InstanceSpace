# Contributing

Thanks for considering a contribution to the Instance Space Analysis (ISA)
Toolkit.

## Setup

See `README.md`'s [Installation Instructions](README.md#installation-instructions)
for the required MATLAB version and toolboxes -- this file doesn't repeat
that, only what changes for development.

No extra setup beyond that is needed: `InstanceSpace.m`/`buildIS.m`/
`exploreIS.m` add `core/`, `output/`, `utils/`, and `deprecated/` to the
MATLAB path automatically the first time any of them runs in a session.

## Before opening a pull request

1. Run `example.m` -- a quick single-configuration smoke test on the
   bundled reference dataset.
2. Run `test_integration.m` -- the exhaustive option-coverage regression
   suite. It prints `EOF:SUCCESS` when every case passes, `EOF:ERROR`
   otherwise (with a per-case summary above it). If your change touches a
   specific option or adds a new one, add or extend a test case rather
   than only relying on the existing ones (see the file's own header
   comment for how cases are structured).

CI (`.github/workflows/tests.yml`) runs both scripts automatically on
every push/PR, but run them locally first too -- catching a failure
before pushing is faster than waiting on a CI run.

## Code style

This codebase has some conventions worth matching rather than
rediscovering:

- Every file starts with a standard header: copyright/SPDX-license block,
  followed by a `Reference:` section citing only the paper(s) that
  specific file actually implements -- not every ISA paper reflexively.
  Copy an existing file's header (e.g. `core/PILOT.m`) as a starting
  point.
- Errors use MATLAB's `identifier` form, namespaced
  `ISA:<FunctionName>:<shortReason>` (e.g. `ISA:ISAvalidateOpts:notInUnitRange`),
  raised with `error(...)` and a message clear enough to act on without
  reading the source -- see `utils/ISAvalidateOpts.m` for the pattern.
- Comments explain *why*, not *what* -- prefer a self-explanatory
  variable/function name over a comment restating the code, and reserve
  comments for a non-obvious constraint, invariant, or the reason behind
  a workaround.
- Options validation is centralised: a new `opts.<section>.<field>` should
  get a default in `utils/ISAdefaults.m` and, where the field has a
  meaningful type/range, a check in `utils/ISAvalidateOpts.m` -- not
  ad-hoc validation inside the algorithm file that consumes it.
- Only `PYTHIA.m`, `TRACE.m`, `PRELIM.m`, and `INIT.m` currently implement
  the train/eval dual-mode convention (an optional trailing
  `trainedModel`/`trainedTrace`/`trainedPrelim` struct argument, dispatched
  on `nargin`); see `CLAUDE.md` and issue #38 before adding this pattern to
  another stage -- `SIFTED`/`PILOT`/`CLOISTER`/`FILTER` deliberately don't
  need it, since none of them recompute anything at explore time.

## Issues and scope

Check the [open issues](https://github.com/andremun/InstanceSpace/issues)
before starting substantial work, in case it's already tracked or
superseded by a larger planned change. For anything beyond a small fix,
opening an issue first (or commenting on an existing one) before writing
code saves rework.
