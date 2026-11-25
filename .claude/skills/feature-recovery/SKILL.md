---
name: feature-recovery
description: >
  Work out what each feature column in a features or metadata CSV actually
  measures, and reproduce it with code that matches the values. Column
  names guide the search; the instance files decide. Use when given a problem
  folder with instances, a CSV of feature columns, and a problem description;
  when asked what a feature means or how it is calculated; or to reproduce a
  feature set whose generating code is missing, undocumented, or written in
  another language. Triggers on "what does this feature measure", "how is this
  feature calculated", "verify these features", "reproduce the feature set",
  "work out the feature definitions", "recover the feature definitions",
  "feature recovery", "feature verification".
---

# Feature recovery

Given a features CSV, the raw instances, and a problem description, work out
what each column measures and write a short Python function that reproduces its
values — a plain-language definition plus code, for every column. Reproduce means
to the precision the column carries, which is often less than it prints.
Descriptive names make this tractable, usually pointing close to the answer, but
a name is a lead, not the thing being established.

## The name narrows the search; the values decide it

Treat a name as the first candidate reading, not as the answer. Accept a
definition because it reproduces the column. When it doesn't, the name has
pointed slightly wrong, and *how* it is wrong says where to look next:

- **Under-specified** — a qualifier is missing. Try weighted, normalised, a
  different denominator, a different log base.
- **Off by a neighbour** — the name's quantity has close relatives: largest for
  smallest, second for first, mean for median, count for fraction.
- **Exhausted** — nothing in the name's family reproduces the column. Stop
  reading the name and search.

## Establish which sources are in scope before you start

Recovery runs in one of two modes, and they are different tasks. Settle which
one you are in before reading anything.

**When the request points at a problem folder**, everything in it is in scope --
and the most valuable thing in it is often not the data. A published feature set
usually has a paper behind it carrying a table that defines every column, and
the folder may hold the extraction code itself. Those settle quickly what a
search over the values settles slowly or not at all: a tolerance, a scaling
constant, a choice between conventions that only a few rows can tell apart. A
constant that no row sits near leaves no fingerprint in the data, and no search
recovers it however long it runs. Find those sources early rather than after a
long search.

**When the request names specific files**, that list is exhaustive. Treat
everything else as absent — including the one file that would settle a question
outright — and say what you would need rather than reading it. A restricted
list usually means the recovery is being used to measure the method rather than
to get the answer, and one extra file ends that with no way to undo it. Widening
the list is the requester's call, never yours; ask.

Either way, note which mode you were in, which files and directories you read,
and where each definition came from. *Output* below says what to do with that.

## Workflow

1. **Align rows to instances, and check both directions**: rows without
   instances, instances without rows. Verification happens on the intersection,
   so report those two counts once and set the unmatched aside — a row with no
   instance and an instance with no row are both unverifiable. Read the counts
   before you do, though: find where the ids actually live rather than assuming
   — the file stem, the parent folder plus stem, or a field inside the file.
   Numbering conventions often differ between subfolders of one dataset, so a
   plausible mapping can bind most rows to the wrong instance while failing
   loudly on only a few, and a large unmatched count is a broken mapping rather
   than missing data.
2. **Profile each column before writing code**: integrality, range, distinct
   count, printed precision, whether values are multiples of `1/n` for some
   instance quantity. A name plus "integer, 1..12, 12 distinct" is a far
   narrower target than the name alone, and the profile constrains the search
   later if the name runs out.

   Sweep the denominator while you are there: for candidate `N`, test whether
   `value * N` is an integer on every row. The `N` that survives names the
   population being counted, which is usually more informative than the count.
   `n` says per element, `n(n-1)/2` says per unordered pair.
3. **Enumerate readings, plural**, from the name, the problem description and
   the profile together. Weighted or not, normalised by what, per-element or
   aggregate, which boundary convention. The enumeration is the work; the first
   guess is rarely the whole answer.
4. **Implement them in a `definitions.py`** exposing `load(path) -> {id:
   instance}` and `DEFINITIONS = {column: callable(instance) -> float}`, then
   `python scripts/verify.py --features features.csv --defs definitions.py`.
   Do the unambiguous columns first: they anchor the instance parsing and the
   alignment, so a mistake there surfaces immediately instead of being blamed on
   a formula. Score every candidate on a small sample first — `--columns` for
   the column under study, `--sample 20` or `30` for the rows — and run the full
   set only once a sample comes back clean. A wrong reading almost always dies in
   the first few rows, so the full run is for confirming a survivor rather than
   for finding the answer.
5. **On a miss, read the verdict**: each one names what kind of wrong the
   reading is, and so which candidate to try next.
6. **Finish when every column is reproduced**, or is one of the few that
   genuinely cannot be — constant columns, and columns that survived the search.

## Two columns that agree on many rows are one measure, twice

Scan every pair of columns for exact equality and count the rows where they
agree. A pair agreeing on many rows and differing on the rest is the same
measure computed twice: once on the instance, once on a transformation of it --
a subset, a rescaled copy, a preprocessed version. The rows where they agree
are the rows where the transformation does nothing, which is what makes the
pattern so recognisable.

This is worth doing early, because it is structural rather than per column: it
can halve the number of distinct measures to recover, and the offsets between
the paired positions tell you how the block is laid out and which measures have
no twin. Identify the transformation on whichever column of the pair constrains
it most tightly.

## Do not reimplement what already exists

Feature sets are usually somebody's published measure family, and those have
maintained implementations. Deriving one by hand from its values is the slow
path, not the opening move.

**Generate candidates wholesale.** Where a package covers the domain's measures,
compute its whole set and match it against every column at once, rather than
deriving one definition at a time. NumPy/SciPy, scikit-learn and networkx cover
most general ground; domain packages exist for most published families and are
worth a search on the measure's name.

**But a library gives you the family, not the conventions.** Defaults rarely
match the original: normalisation, missing values, multiclass decomposition,
boundary and tie rules all still have to be pinned by verification. It is still
a large win — it collapses the search from "what is this quantity" to "which
options were set".

**And the set itself may be published**, in which case the paper defining it
beats any search — see *Establish which sources are in scope* above for whether
that one is yours to read.

Bring dependencies in freely for the search; keep the shipped `definitions.py`
lean where a few lines do the job.

## Keep the search loop cheap

What costs time here is not one definition running slowly. It is the loop —
readings enumerated, most of them discarded, each scored against every instance.
Speeding up a candidate you are about to throw away is the wrong end of it, and
rewriting one in a faster language is the wrong end twice over. Cut the loop
first; only then look at the survivor.

**Score a candidate on a small sample first, then confirm on everything.**
`--sample 30` scores thirty rows, `--columns` restricts to the column under
study. Twenty or thirty rows is the working size: a structurally wrong reading
usually dies on the first two or three, so nearly every candidate can be killed
without touching the corpus. Once a sample comes back clean, rerun that column on
the full set — a sample refutes a reading and never confirms one.

The sample is drawn deliberately, not at random, because what it has to catch is
a missing branch. Those sit at the ends of a column's value range — the
degenerate instance, the row sitting exactly on a threshold — so the sample takes
each column's global extremes first. If the CSV records a source for each row, it
then spreads the remaining picks across those sources rather than letting them
cluster in one, since the other common missing branch is a whole generator
handled differently. The run prints what the sample covered.

What a sample cannot do is find a *rare* error, and no sample size worth running
fixes that: as the disagreement narrows towards a single row, the sample you
would need approaches the whole corpus. So `--sample` is for killing wrong
readings, which is most of the search. Characterising a handful of failing rows —
the work under *Reading a PARTIAL* — is full-run work, and reading a clean sample
as though it were a clean run is how a missing branch survives to the write-up.

Two cautions. Refutation survives a tiny sample but the verdict *label* does not
— a real miss can report as `MONOTONE` rather than `PARTIAL` — so read a small
sample as pass/fail and nothing finer. And a column can be flat within a sample
while varying over the CSV; that is reported as `FLAT`, not `CONSTANT`, since the
two mean opposite things.

**Compute shared intermediates once.** `verify.py` calls one definition per
column per instance, so anything several columns need — a pairwise distance
matrix, a sorted array, an expanded element list — is recomputed once per column
unless the instance object holds it. Cache on the instance, not inside the
definition:

```python
class Inst:
    @cached_property
    def elements(self):       # a dozen columns below share this
        ...
```

`load` may also take an optional `ids` keyword; when `--sample` is in use it
receives the ids actually being scored, so only those instances are parsed. Worth
having when parsing is the expensive half, which it often is.

**Measure before optimising either.** The run prints a timing line naming the
slowest columns and the per-instance cost. A guess about where the time went is
usually wrong, and the answer decides whether any of this is worth doing.

## Verdicts

| Verdict | Meaning | What to do |
| --- | --- | --- |
| `EXACT` | reproduces every aligned row | done |
| `EXACT~` | agrees within float noise | done; iterative and spectral methods cannot reproduce bit-for-bit |
| `EXACT~f32` | every row agrees to ~7 significant digits | done; the column was stored at single precision, or accumulated long |
| `AFFINE` | right quantity, wrong scale | the reported slope names the missing step — `1/n`, `1/ln 2`, `2/(n(n-1))` |
| `MONOTONE` | right ordering, unknown transform | try transform variants only, not new formulas |
| `PARTIAL` | matches some rows | characterise the failures, then refine this reading |
| `NO MATCH` | this reading is unrelated to the column | discard it and try the next candidate; if the name is exhausted, search |
| `CONSTANT` | column never varies | unidentifiable from values; record it and move on |
| `FLAT` | constant within a `--sample`, but not in the column | not a verdict; raise `--sample` or drop it |

The three `EXACT` verdicts are answers. `AFFINE` and `MONOTONE` are not: a
definition that merely correlates says nothing about how the column was
calculated, however well the name fits.

Both are reserved for a reading that misses on most rows. One that is exact on
most of the non-zero rows and far out on the rest reports as `PARTIAL`, however
close the correlation, because what it lacks is a branch, not a transform.

They are leads rather than discards, though. Both say the quantity is probably
right and only the transform is wrong, which is a far narrower search than
starting over — `NO MATCH` is the one verdict that means throw the reading away.
On a small sample neither label is trustworthy, so rerun the column on the full
set before acting on it: a missing branch that happens to preserve ordering
across the sampled rows reports as `MONOTONE` and resolves to `PARTIAL` once the
full set is scored. Which one it settles on decides the next move — transform
variants for a `MONOTONE` that survives, failure characterisation for a
`PARTIAL`, and for a bounded column the squash under *A bounded column may be a
squashed unbounded one* is the transform to try first.

## Rounding is not a miss

A CSV usually prints more decimals than the number behind it ever had. A column
computed in single precision and printed to nine decimals looks precise and is
not — `0.699999988` is `float32(0.7)`, and no correct reimplementation in double
precision will land closer than about `1e-7` relative. Demanding the printed
decimals there rejects the right answer.

These columns are downstream inputs to an instance space, where the projection is
computed from standardised features: a deviation in the seventh significant digit
cannot move an instance visibly. A *wrong definition* does, which is why `AFFINE`
and `MONOTONE` are never answers, however close the correlation.

So on a near miss the question is not how large the error is but **whether it is
the same size on every row**:

- Small everywhere → rounding. Take `EXACT~f32` and move on.
- Small on most rows, far out on a few → a missing branch, and the far-out rows
  are the ones to characterise. `verify.py` reports both counts and diagnoses the
  failures at the looser bound, so a handful of genuine misses is not buried
  under hundreds of rounding-level ones.

Storage width cannot be recovered once a CSV has rounded to fixed decimal places,
so do not try to detect it — score at both bounds and report which one the column
meets. That is a complete result, not a hedged one.

## Reading a PARTIAL

**Scattered failures are edge cases; a predictable set is a missing branch.**
Both look identical as a count. Three stray rows are worth logging and shipping;
forty that all share a property mean the implementation lacks a branch, however
good the overall match rate looks. Test whether failure is predictable from
instance properties, and weigh any pattern against its base rate — "every
failure is a dense instance" means nothing when most instances are dense.

**Boundary conventions are the usual culprit.** Strict versus non-strict differs
only on values sitting exactly on the threshold: few enough rows to look like
noise, and enough of them to be the wrong definition. Sweep the conventions
rather than assuming one, and do not assume classes that look like a partition
actually are — they may overlap, or leave a gap.

**A contradiction refutes the reading entirely.** If the same input is counted
one way in one row and differently in another, no deterministic function of that
input reproduces the column: either something outside the instance files is
involved, or the original computed an intermediate differently — an extra
normalisation step changes rounding without changing the intended quantity.

**Some columns are not reproducible, and saying so is a result.** When the
failures are confined to instances sharing a structural property and no
definitional variant fits, the difference may lie in the original's
implementation rather than in its definition. Two checks make that a conclusion
rather than a shrug. First, confirm the property is *necessary*: if every
failing row has it and no row without it fails, its handling is the whole
story. Second, find out whether the published value depends on the instance's
values alone, by comparing rows whose instances hold the same values stored
differently. If they agree, the column is computable from the data and the
search is not finished — go back to the failing rows. If they disagree, the
original distinguishes elements that are indistinguishable by value, and
nothing derived from values can match it. If the corpus has no such rows, the
check cannot be run; say so rather than count it as passed.

Only when both checks point at the implementation, report the column open with
the list of readings ruled out, and stop. Fitting the column to reproduce the
corpus would give code that matches these rows and computes nothing on the next
instance, which is the opposite of the point.

## A bounded column may be a squashed unbounded one

Feature sets destined for an instance space are often bounded into `[0, 1]`
before they are stored, because the pipeline standardises every column. When the
natural quantity is unbounded and the column is not, look for the squash and
invert it before searching. `arctan(x/D)/(pi/2)`, `tanh(x/D)`, `x/(x+D)`, a
logistic and a plain cap are all in use; try each rather than assuming one.

The tells are in the profile. A column that approaches its bound without
reaching it, a maximum sitting exactly on a round value, or a bounded column
whose name promises a ratio, a spread or a gap — all unbounded quantities.
Invert each candidate squash and look at what comes back: raw values are far
easier to recognise than squashed ones, and they often turn out to be integers
or simple ratios the squash had hidden.

A squash usually carries a scale `D`, a round constant chosen so the mapping is
roughly linear over the range that occurs, and it will not fall out
of the values on its own. Recover it by solving: take the inverted value, divide
by the quantity you think it measures, and check whether the result is constant
across the corpus. A constant is the answer; a divisor that varies means the
reading is wrong.

A ratio can be written either way round — `a/b` or `b/a` — and the column's
name will not always tell you which. Check the direction against the values
before fitting, or every candidate formula is off by an inversion.

## Searching, once the name is exhausted

A dead name means the column is unnamed, not undefinable. Compute a broad set of
cheap instance quantities — moments of every natural array, the same normalised
by instance size, ratios between them, counts past canonical thresholds — and
look for one that reproduces the column exactly. Sweep a library's whole measure
set here too, in both weighted and unweighted form: instance files often carry
weights the problem statement never mentions. Use the column profile to prune —
an integer column bounded by `n` rules out most of the space before anything is
computed.

**Ask what the column is a function of, before asking what it equals.** Pick a
candidate quantity, group the rows by it, and check whether the column is
constant within every group. This is cheap, needs no formula, and is decisive in
both directions: if one integer determines the column on every group, the search
has collapsed from all of feature space to one variable, and a closed form is
usually visible in a table of the two. If nothing determines it, whole families
are ruled out at once. Run the test before enumerating readings, not after.

A column that is a function of one integer is also a column carrying no more
information than that integer — worth saying in the write-up, because a feature
set that encodes the same quantity twice will distort the projection.

## Output

One self-contained function per column, standard library plus NumPy/SciPy, the
plain-language definition in its docstring, and the verification result recorded
beside it — with both the denominator, which is the aligned subset, and the
precision it holds at: "exact on *k* of the *m* rows that have instances", or
"exact at single precision on *k* of *m*".

Score each column at both bounds, as *Rounding is not a miss* says — but count
the disagreeing *rows* at one fixed bound across the whole set. Counting each row
at its own printed precision instead flatters the columns that print few
decimals — a row printed as `0` tolerates an error of half a unit — so two
columns quoted that way are not comparable, and a column can look several times
better than it is.

Say how well each definition is pinned down, because the next person needs to
know which parts they can build on. Three levels are worth telling apart: what
the values reproduce, which the verification result already records; what a
source supplied and the values then confirmed; and what nothing in the data
settles — a tolerance inside a threshold, a tie rule, an imputation choice —
taken on a source's word or chosen by you. That last level is the one that has
to be labelled, because written as code it looks exactly like the other two.

Then the residue: columns nothing reproduced, columns that are constant, and —
as a note for whoever relies on the names later — columns whose name turned out
to be a poor description of what they measure.
