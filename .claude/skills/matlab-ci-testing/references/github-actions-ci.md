# GitHub Actions CI for MATLAB

## Minimal working workflow

```yaml
name: Tests

on:
  push: {}
  pull_request:
    branches: [master]
  workflow_dispatch: {}

jobs:
  test:
    runs-on: ubuntu-latest
    # Observed run: setup-matlab's system-dependency apt-get step stalled on
    # a runner-side network/mirror hang for ~2h45m with no output at all
    # (unrelated to any code here) before it had to be cancelled manually.
    # Caps that failure mode to a fast, visible failure instead of silently
    # running toward GitHub's 6-hour default job timeout.
    timeout-minutes: 90
    steps:
      - uses: actions/checkout@v4

      - uses: matlab-actions/setup-matlab@v2
        with:
          # Pinned, not 'latest': the README/docs promise support from a
          # specific minimum MATLAB release onward, so CI should keep
          # validating that minimum specifically rather than silently
          # drifting to whatever release 'latest' resolves to as new
          # MATLAB versions ship. A regression only reachable on an older
          # release would otherwise go undetected once 'latest' moves past
          # it.
          release: R2025a
          products: >
            Global_Optimization_Toolbox
            Optimization_Toolbox
            Parallel_Computing_Toolbox
            Statistics_and_Machine_Learning_Toolbox
            Financial_Toolbox

      # Both scripts below abort on any uncaught error, which run-command
      # turns into a nonzero exit code and a failed step -- no separate
      # pass/fail parsing needed here.
      - uses: matlab-actions/run-command@v2
        with:
          command: example, test_integration

      - uses: actions/upload-artifact@v4
        if: failure()
        with:
          name: test-outputs
          path: test/data/
          retention-days: 14

      # Coverage report is small and useful either way, not just on
      # failure.
      - uses: actions/upload-artifact@v4
        if: always()
        with:
          name: coverage-report
          path: coverage.xml
          retention-days: 14
```

## Discovering the right toolbox list

There's no shortcut for this beyond running CI and reading what breaks.
Start with an empty or minimal `products:` list, run CI, and read the
first error:

- `Undefined function 'X' for input arguments of type ...` where `X` is
  a toolbox function (not your own code) → add the toolbox that ships
  `X`. In this repo, `PRELIM.m`'s auto-normalisation calls `boxcox()`,
  which isn't obvious from the function name alone that it needs the
  **Financial Toolbox** specifically (not Statistics and Machine
  Learning, where a naive guess would land) — confirmed only by making
  CI fail without it, then adding it and watching the same failure
  disappear.
- A licensing/entitlement error naming a toolbox → same fix, different
  error shape.

Repeat until a full run gets past `setup-matlab` and into your actual
test command. Don't guess a "complete" list up front — it's faster to
let CI tell you exactly what's missing, one iteration at a time, than
to over-provision products you don't actually need.

## Why pin a release instead of `latest`

`release: latest` seems convenient, but it means CI validates whatever
MATLAB happens to be current *today* — not the minimum version your
documentation promises to support. If a user on your documented minimum
release hits a bug that only reproduces on an older MATLAB (a function
signature change, a default-value change, a toolbox behavior change),
`latest` will never catch it, and it'll silently drift further from your
actual support matrix every time a new MATLAB version ships. Pin to the
specific release named in your README/install docs instead — if you
later raise the minimum supported version, update the pin deliberately,
as a decision, not as a side effect of time passing. `matlab-actions/setup-matlab@v2`
accepts a specific release string (e.g. `R2025a`) exactly like `latest`,
so this is a one-line change with no other tradeoffs.

## The timeout is a safety net for a specific observed failure mode, not decoration

GitHub's default job timeout is 6 hours. A runner-side network/mirror
hang during `setup-matlab`'s system-dependency installation step can
stall with **zero log output** for hours before anything fails —
indistinguishable from "still working" unless someone happens to notice
and cancel it manually. A much shorter `timeout-minutes` (90 was ample
margin over this repo's actual ~30-45 minute successful run time) turns
that failure mode into something that fails fast and visibly instead of
silently burning CI minutes — and, more importantly, instead of leaving
a human waiting on a run that will never finish on its own. Set this
based on your own observed successful run time plus real margin, not a
guess.

## Reading a failing run's logs

`matlab.unittest`'s runner prints two things worth knowing about at the
*end* of a run (both survive being near the tail even in a huge log):

- A **Failure Summary** table (`Name | Failed | Incomplete | Reason(s)`)
  listing every failing test by name.
- A final `[TEST] [PASS]/[FAIL] <TestClass>/<testMethod>` line per case
  if your runner prints one (worth adding — see
  `matlab-unittest-architecture.md`).

If a log-retrieval tool caps how much content it returns (common: it
returns roughly the last N kilobytes of the log regardless of how much
you ask for), and the failing test ran *early* in a long suite with lots
of verbose output after it, the actual stack trace/error detail for that
test may simply be unreachable through that tool. When that happens,
don't keep re-requesting larger windows hoping to get lucky — reason
from what you *can* see (the Failure Summary's test name, the PASS/FAIL
line) plus direct source-code reading of that specific test and the code
it exercises. In practice this is usually enough to pinpoint the bug
without the literal exception text.

One more distinction worth knowing when reasoning about a failure:
`matlab.unittest` reports each test as **Errored** or **Failed**, and
they mean different things. *Failed* means a `verifyXxx` qualification's
condition was false — an assertion mismatch. *Errored* means MATLAB
threw an uncaught exception — a real crash, not a value comparison. If a
log shows "Errored" with no other detail reachable, look for something
that would throw outright: e.g. a test referencing a class property that
turns out to be `Access = private` (from *outside* the class, this
throws an access-denied `MException`, not a comparison failure) is a
concrete example that produced exactly this symptom in this repo.

## Confirm the run you expect actually exists

A `git push` can occasionally hit a ref-lock race (typically visible as
`cannot lock ref '...': is at <sha> but expected <sha>`) where the push
still lands the commit on the remote successfully, but the corresponding
push-triggered workflow run never fires. `git log`/`git fetch` showing
the commit present on the remote is *not* sufficient confirmation that
CI ran for it — separately check that a workflow run exists with that
exact commit SHA as its `head_sha` before treating silence as "still
running" or, worse, assuming the previous green run still covers it. If
no run exists, the fix is usually just to push again (a trivial
follow-up commit, or the next real commit, will trigger normally).
