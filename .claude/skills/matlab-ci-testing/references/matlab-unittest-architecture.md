# matlab.unittest architecture

## From a bare script to TestCase classes

**Before**: one script, a cell array of hand-built case structs, a
`try`/`catch` around each, manual pass/fail bookkeeping printed at the
end.

**After**: a `tests/` folder of `matlab.unittest.TestCase` subclasses,
one per concern (not one giant class) — e.g. option-surface coverage in
one class, staged class-API round-trips in another, a legacy-migration
table in a third, targeted bug regressions in a fourth. The original
script becomes a thin runner:

```matlab
import matlab.unittest.TestSuite
import matlab.unittest.TestRunner
import matlab.unittest.plugins.CodeCoveragePlugin
import matlab.unittest.plugins.codecoverage.CoberturaFormat

repoRoot = fileparts(mfilename('fullpath'));   % see "Path resolution" below
addpath(repoRoot);

testDataDir = fullfile(repoRoot, 'test', 'data');
if ~isfolder(testDataDir)
    mkdir(testDataDir);
end

suite = TestSuite.fromFolder(fullfile(repoRoot, 'tests'), 'IncludingSubfolders', true);
runner = TestRunner.withTextOutput();

coverageReportFile = fullfile(repoRoot, 'coverage.xml');
sourceFolders = {repoRoot, fullfile(repoRoot, 'core'), fullfile(repoRoot, 'output'), fullfile(repoRoot, 'utils')};
runner.addPlugin(CodeCoveragePlugin.forFolder(sourceFolders, ...
    'IncludingSubfolders', false, 'Producing', CoberturaFormat(coverageReportFile)));

results = runner.run(suite);

fprintf('\n[TEST] ================= Summary =================\n');
for i = 1:numel(results)
    if results(i).Passed, status = 'PASS'; else, status = 'FAIL'; end
    fprintf('[TEST] [%s] %s\n', status, results(i).Name);
end
nPassed = sum([results.Passed]);
nCases  = numel(results);
fprintf('[TEST] %d/%d cases passed.\n', nPassed, nCases);

if nPassed == nCases
    fprintf('EOF:SUCCESS\n');
else
    fprintf('EOF:ERROR\n');
    error('MyProj:testRunner:caseFailures', '%d of %d test cases failed.', nCases-nPassed, nCases);
end
```

Two things in that last block matter beyond style:

- **Explicitly `error(...)` on failure.** `TestRunner.run` doesn't itself
  make the calling script exit nonzero — printing "FAIL" to the console
  and returning normally would leave `matlab-actions/run-command`
  thinking the step succeeded. Raising an error is what actually turns a
  test failure into a red CI check.
- **Keep any pre-existing sentinel convention** (here, `EOF:SUCCESS`/
  `EOF:ERROR`) if other scripts or external tooling already scrape for
  it. Check for that *before* changing how the runner reports its
  result — it's easy to "clean up" old bookkeeping that's still a load-
  bearing contract for something outside the file you're editing.

## TestSuite discovery, TestParameter, and CodeCoveragePlugin

- `TestSuite.fromFolder(folder, 'IncludingSubfolders', true)`
  auto-discovers every `TestCase` subclass under `folder`. No manual
  registration list to keep in sync as tests are added.
- `TestRunner.withTextOutput()` streams readable progress to the
  console as tests run — useful for a CI log a human might actually read
  while a run is in progress, not just after it finishes.
- Parameterize instead of hand-rolling a cell array of cases:

  ```matlab
  properties (TestParameter)
      OptionCase = pipelineOptionCases();  % local function returning a struct-of-structs
  end

  methods (Test)
      function testOption(testCase, OptionCase)
          % runs once per field of the struct OptionCase returns
      end
  end
  ```

  Each field becomes its own test instance, reported as
  `ClassName/testOption(OptionCase=fieldName)` — you get per-case
  pass/fail and the standard diagnostic for free, instead of manually
  tracking which case failed. `pipelineOptionCases()` can live as a
  plain local function at the bottom of the same file (MATLAB allows
  local functions after a `classdef` block).
- `CodeCoveragePlugin.forFolder(sourceFolders, 'IncludingSubfolders', false, 'Producing', CoberturaFormat(file))`
  produces a standard Cobertura XML report — consumable by most CI
  coverage dashboards/badges without extra tooling. Pass every source
  folder you actually want measured; it won't infer them from what the
  tests happen to `addpath`.

## A non-obvious TestClassSetup gotcha

`methods (TestClassSetup)` runs once per class, not once per test
method — but the important, easy-to-get-wrong part is that **the same
class instance is reused across every `(Test)` method in that class**.
A property set in `TestClassSetup` (a shared fixture directory, a
built-once model) is visible in every test method afterward, and if one
test method *mutates* an instance property, that mutation is visible to
whichever test method happens to run next. This is not the isolated,
fresh-instance-per-test model some other frameworks default to.

Two consequences worth designing around:

- Don't assume execution order is something you can rely on formally —
  it isn't documented/guaranteed — but tests observably tend to run in
  file-declaration order within a class. If several checks genuinely
  depend on a specific sequence (e.g. a staged build-then-explore
  round trip), keep that sequence as *one* test method rather than
  splitting it across several methods hoping they run in the order you
  wrote them.
- If two independent test methods both need a fresh, uncorrupted
  fixture (e.g. both build their own model from the same shared case
  directory), build each one from scratch inside its own method rather
  than relying on state another method might have already mutated.

## Path resolution: the working directory is not guaranteed

This is the single most disruptive difference between "run this script
directly" and "run these same files through `matlab.unittest`": **the
working directory during test execution is not guaranteed to be wherever
you launched the runner from.** Observed directly in CI: `TestClassSetup`
methods ran with the current folder set to `tests/` itself, turning a
relative literal `'./test/data/...'` into `'tests/test/data/...'` and
failing every `copyfile`/`readtable` call with a plain "file not found."

The same root cause is broader than test-data paths. Any file that's
only resolvable via MATLAB's *implicit current-folder lookup* — never
formally `addpath`-ed, just "happens to work" because a normal script's
cwd was the repo root — breaks the same way, surfacing as "Undefined
function" instead of "file not found." A top-level entry-point script
sitting at the repo root, relying on the fact that scripts are usually
launched from there, is a common instance of this.

### The fix: resolve paths from the file's own location, not the caller's cwd

`mfilename('fullpath')` always returns the calling file's own physical
location on disk, regardless of what the current working directory
happens to be. Build a small helper around it:

```matlab
function root = testRepoRoot()
% One level up per directory this file sits below the repo root.
root = fileparts(fileparts(mfilename('fullpath')));   % e.g. tests/testRepoRoot.m -> repo root
if ~(endsWith(root, '/') || endsWith(root, '\'))
    root = [root filesep];
end
end
```

A file sitting directly at the repo root (like the thin runner script
above) needs only one `fileparts()` call:
`repoRoot = fileparts(mfilename('fullpath'));` — count the actual
directory depth, don't copy-paste the double-`fileparts` version
unthinkingly.

### Apply it everywhere, not just at the first failure

This is the part that's easy to under-scope: fixing *one* relative path
does not fix the class of bug. Every independent relative-path-shaped
call site breaks on its own, and in this session each one only surfaced
as its own separate CI failure, one push at a time — three rounds of
"fix the one that just failed, push, watch the *next* one fail." Grep
for every relative-path-looking literal and path-taking-function call in
one pass instead:

- `addpath(...)` — for any repo-root-level file not already formally on
  the path (`addpath(repoRoot)`, not `addpath(pwd)` or omitted).
- Test-data / fixture directories — `fullfile(repoRoot, 'test', 'data')`,
  not `'./test/data/'`. Each test *class* that builds its own fixture
  directory needs this independently (e.g. a `Constant` property
  `RootDir = [testRepoRoot() 'test/data/'];`), not just the top-level
  runner.
- `TestSuite.fromFolder(...)` — pass `fullfile(repoRoot, 'tests')`, not
  the bare string `'tests'`.
- `CodeCoveragePlugin.forFolder(sourceFolders, ...)` — build
  `sourceFolders` from `repoRoot`, not relative literals like
  `{'.', 'core', 'output'}`.
- The coverage report path itself — `fullfile(repoRoot, 'coverage.xml')`,
  so it lands somewhere predictable regardless of cwd.
- `mkdir(...)` for any directory the test setup creates on demand.

A reasonable way to sanity-check you got all of them: search the whole
`tests/` folder (and the runner script) for both `'./'`-prefixed string
literals and any call to a path-taking function (`addpath`, `mkdir`,
`fromFolder`, `forFolder`, `isfolder`, `copyfile`, `readtable`, ...)
passed something that isn't already an absolute/`fullfile`-built path.
