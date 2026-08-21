# MATLAB test-writing patterns

## `assumeTrue`/`assumeFalse` as precondition guards, not assertions

`matlab.unittest` qualifications come in families with different
severities: `verifyXxx` (soft — records a failure, keeps running the
rest of the method), `assertXxx` (hard — stops the method immediately on
failure), and `assumeXxx` (a *precondition* — if false, the test is
marked **skipped**, not failed or passed, and stops there).

Use `assumeTrue`/`assumeFalse` when a test needs a specific precondition
to even be meaningful, and that precondition isn't something the test
itself controls precisely — e.g. confirming a randomly-constructed
`alphaShape` actually produced multiple disconnected regions before
testing multi-region behavior, or confirming a model build actually
produced at least two trained algorithms before testing a code path that
needs that. Without this, a test whose setup happens not to hit the
scenario it's meant to exercise either fails misleadingly (looks like a
real regression when it's actually an unlucky data draw) or — worse —
silently no-ops and reports a false pass. `assumeTrue` reports the
correct third outcome: "inconclusive, not this run's fault."

Prefer `verifyXxx` for the actual assertions about behavior within one
test method, since a single method often checks several related
properties and you want to see *all* of them reported, not just the
first failure. Reach for `verifyWarning(@() fcn(...), 'My:warning:id',
msg)` specifically when checking that a **warning** (not an error) fired
with a given identifier — MATLAB warnings don't throw, so `verifyError`
doesn't apply to them.

## Testing an internal/helper function via workspace injection

If a codebase already shares helper functions across several top-level
scripts by injecting them into the caller's workspace (a common MATLAB
pattern: a script-style file defines local functions, then calls
`assignin('caller', 'name', @name)` for each one it wants to expose), a
test method can use the exact same mechanism to reach functions that
would otherwise be untestable — not on the path, not exposed any other
way:

```matlab
function testSomeHelper(testCase)
    scriptfcn;   % injects everything scriptfcn.m exposes into this method's workspace
    result = someInjectedHelper(args);   % callable as if it were a local function
    testCase.verifyEqual(result, expected);
end
```

If the function you actually want to unit-test is one level *below* what
the file currently exposes (e.g. it exposes a top-level
`traceAlphaBoundary`, which internally calls a lower-level
`traceOneRegion` you want to test directly, with tighter control over
the exact input than going through the top-level function would allow),
just add another `assignin('caller', 'traceOneRegion', @traceOneRegion);`
line. It's a one-line addition, and it's a far better lever for testing
a specific inner behavior than reconstructing whatever complex
input/geometry the outer function normally needs to reach it indirectly.

## Prefer synthetic, hand-built edge cases over fragile geometry/randomness

When the function under test operates on a well-defined intermediate
data structure — not a black-box object from some upstream algorithm —
construct that structure by hand for the exact edge case you want,
rather than trying to coax a "realistic" data-generating process into
producing it.

Concrete example: testing how a boundary tracer handles "a region with a
hole" (two disconnected boundary cycles) could be attempted by building
a real `alphaShape` shaped like an annulus and hoping MATLAB's own
region-detection algorithm produces the intended internal structure —
fragile, since it depends on undocumented `alphaShape` internals, needs
tuning of alpha/point-spacing to reliably reproduce, and could silently
break on a future MATLAB version. The function under test actually
operates on a simple edge-list graph (pairs of vertex indices plus
coordinates), so construct that graph directly instead:

```matlab
bf = [1 2; 2 3; 3 4; 4 1; ...   % outer 4-cycle
      5 6; 6 7; 7 5];           % disconnected inner 3-cycle (the "hole")
bv = [0 0; 4 0; 4 4; 0 4; 1.5 1.5; 2.5 1.5; 2 2.5];
verts = traceOneRegion(bf, bv);
```

This gives full, deterministic control over exactly the scenario being
tested, independent of how some upstream process happens to produce that
topology in the wild. General principle: test at the level of the
simplest well-defined representation the function actually consumes, not
only at the "realistic" entry point real callers use — and expose that
level for testing (via the injection pattern above, if needed) when it
isn't already reachable.

When randomness genuinely is what's under test (e.g. confirming a seed
makes two runs reproducible), seed it explicitly at the top of the test
(`rng(1, 'twister');`) so the *test itself* stays deterministic even
though it's exercising the system under test's own random behavior.

## Callback testing: the handle-class accumulator trick

MATLAB anonymous functions (`@(args) expr`) can only be a single
*expression* — they cannot contain assignment statements. So a callback
like `@(stageName, result) log{end+1} = stageName` is a syntax error,
not valid MATLAB, which makes "record what a callback was invoked with"
look harder than it is.

The fix: `matlab.unittest.TestCase` is itself a **handle class** (mutable
by reference). Give the test class its own accumulator properties and a
plain instance method that appends to them:

```matlab
properties
    StageLog = {}
end

methods
    function logStage(testCase, stageName)
        testCase.StageLog{end+1} = stageName;   % a statement -- fine in a real method
    end
end

methods (Test)
    function testCallbackFires(testCase)
        testCase.StageLog = {};
        cb = @(stageName) testCase.logStage(stageName);   % single expression -- valid anonymous fn
        systemUnderTest(cb);
        testCase.verifyEqual(testCase.StageLog, {'expectedStage1', 'expectedStage2'});
    end
end
```

The anonymous function's body is a single expression — a method call —
even though that method call has side effects (mutating the handle
object) inside it. Reset the accumulator property at the start of each
scenario within a test method if you're exercising the same callback
against more than one code path and expect different accumulated results
for each.

## Name and comment regression tests as documentation of the bug

Name each regression test for the scenario/guarantee it protects, not
generically — `testSeedReproducibility`, `testTraceAlphaBoundaryMultiRegion`,
`testPrelimEvalModeAlgoAlignmentAfterPruning`, not `testBugFix3`. In the
leading comment, state:

1. What issue/review round it came from, if applicable.
2. The actual root cause — not "this used to be broken," but *why* it
   was broken (what assumption failed, under what condition).
3. Why the test is constructed the specific way it is, especially if
   that construction looks like an arbitrary choice at a glance — e.g.
   "prunes the *middle* algorithm, not the last one, since that's what
   exposes the positional-misalignment half of the bug" or "shortens
   *both* arrays by the same amount, since a mismatch-count check alone
   wouldn't catch this." A future reader (human or another agent)
   shouldn't have to reverse-engineer why the test looks oddly specific.

This has a second, load-bearing purpose beyond documentation: it makes
it possible to tell, before writing a new test, whether a specific
historical failure mode is already covered.

## Validate a new regression test against the *old*, buggy code

Before trusting that a new regression test actually protects against
what it claims to, check that it would have **failed** on the code as it
was before the fix — not just that it passes now. It's easy to write a
test that technically exercises the fixed code path but doesn't actually
discriminate between the buggy and fixed behavior. A concrete instance
from real use: a fix added a warning for a previously-silent failure
mode; the *first* version of that fix didn't actually change the
condition the warning was meant to catch, so the warning still never
fired — and a test that only checked "did a warning fire" (without also
checking the returned value's actual content) would have passed either
way, giving false confidence that the fix worked. Checking the test
against the pre-fix code — by re-reading the old logic and reasoning
through it, or by literally reverting the fix locally and re-running —
catches this before it ships as a green checkmark over a bug that's
still there.
