---
name: matlab-ci-testing
description: >
  Set up or debug a MATLAB continuous-integration pipeline and
  matlab.unittest test-suite architecture on GitHub Actions --
  matlab-actions/setup-matlab + run-command configuration (toolbox
  dependency discovery, pinning a MATLAB release vs 'latest', job
  timeouts as a safety net against runner hangs), migrating a hand-rolled
  test script to matlab.unittest.TestCase classes (TestSuite.fromFolder,
  TestParameter, CodeCoveragePlugin/CoberturaFormat), and the classic
  repo-root path-resolution bug where matlab.unittest silently changes
  the working directory during test execution. Use this whenever the
  user is adding CI to a MATLAB repository, migrating MATLAB tests to
  matlab.unittest, debugging a failing MATLAB GitHub Actions run, or
  hitting "file not found"/"undefined function" errors that only happen
  inside matlab.unittest but not when running a script directly. Also
  covers test-writing patterns specific to MATLAB (assumeTrue
  precondition guards, testing internal helpers via assignin injection,
  synthetic edge-case construction over fragile geometry, callback
  testing via handle-class TestCase properties) and how to read a
  matlab.unittest CI failure log efficiently.
---

# MATLAB CI + matlab.unittest architecture

Distilled from building this out for real on `andremun/InstanceSpace`:
adding GitHub Actions CI to a MATLAB repo that had none, migrating its
hand-rolled test script to `matlab.unittest`, and then living through
five rounds of code review that each found real, previously-invisible
bugs — several of them *in the test infrastructure itself*, not just the
code it was testing. That last part is not a footnote: a test framework
you can't trust is worse than no tests, because it actively hides bugs
behind a green checkmark. The patterns here exist because a plausible
first attempt at each of them was wrong in a specific, discoverable way.

## When this applies

- The repo has no CI, or CI exists but doesn't run MATLAB tests.
- Tests are a bare script (`try`/`catch` per case, hand-rolled
  pass/fail bookkeeping) and need to move to `matlab.unittest`.
- A MATLAB GitHub Actions run is failing and the cause isn't obvious.
- Code that works fine run directly in MATLAB throws "Undefined
  function" or "file not found" only when run through `matlab.unittest`.
- Writing a MATLAB regression test for a specific bug, especially one
  involving randomness, an internal/private helper, or a callback.

## The four pieces, and where to read about each

1. **CI workflow itself** — `matlab-actions/setup-matlab`,
   `matlab-actions/run-command`, which MATLAB toolboxes to declare, why
   to pin a specific release instead of `latest`, and a timeout as a
   safety net. See `references/github-actions-ci.md` — it also covers
   how to actually read a failing run's logs without fighting a
   log-retrieval tool's size cap, and a real git-push race that silently
   dropped a CI run for one commit.

2. **matlab.unittest architecture** — replacing a bare script with
   `TestCase` classes, `TestSuite.fromFolder`, parameterized tests via
   `TestParameter`, and `CodeCoveragePlugin`. See
   `references/matlab-unittest-architecture.md`.

3. **The repo-root path problem** — `matlab.unittest` does not guarantee
   the working directory stays where you launched it from, so any
   relative path silently breaks, one call site at a time. This is
   folded into `references/matlab-unittest-architecture.md` too (it's a
   consequence of the architecture, not a separate topic), with the
   fix pattern and — importantly — the checklist of *every* place it
   needs applying, not just the first one you find.

4. **Test-writing patterns** — `assumeTrue` vs `verify*`, testing an
   internal helper by injecting it into the test's workspace, building
   a synthetic edge case by hand instead of fighting geometry/randomness
   into producing one, the handle-class accumulator trick for testing
   callbacks (MATLAB anonymous functions can't contain assignments), and
   how to name/comment a regression test so it documents the bug it
   guards against. See `references/test-writing-patterns.md`.

## The one meta-lesson worth internalizing

Every non-trivial fix in this area — the path-resolution helper, the
coverage plugin wiring, a regression test itself — got something wrong
on the first pass, and the mistake was only caught by either a
subsequent CI failure or a second, more careful look. Two habits from
this session are worth carrying forward on any MATLAB CI/test work:

- **Fix the class of bug, not the instance.** When a relative path
  breaks inside `matlab.unittest`, don't fix that one call site and
  move on — grep for every other relative-path-shaped literal in the
  same file (and sibling files following the same pattern) and fix them
  together. The debugging cycle in this session found three separate
  relative-path breakages across three separate CI runs, each one
  discoverable in the first pass with a bit more suspicion.
- **A regression test can itself be wrong.** Before trusting a new
  regression test, check whether it would actually have failed against
  the *old*, buggy code — not just that it passes against the fix. One
  fix in this session added a warning for a previously-silent failure
  mode; the first version of that fix didn't actually change the
  code path the warning was meant to catch, and the warning still never
  fired. A test that only checks "did a warning fire" without also
  checking the returned *content* would have passed anyway.
