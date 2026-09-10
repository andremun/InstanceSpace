# Documentation site build status

Tracks issue #53's 29-page inventory against what's actually built. See
`llm-server/docs/instancespace_docs_site_framework.md` for the full
framework spec this build follows (task registry, categories, gates).

Toolchain (§4a, confirmed 2026-09-09): Markdown-first source
(`doc/src/*.md`) generated to HTML (`doc/html/*.html`) via
`doc/generate.py` (python-markdown), sharing `doc/style.css`. Nested
`opts.*`/`out.*` struct fields use dotted-key flattening (own row,
lightly indented), not a nested table.

**Everything below was built through the real Planner (`gemma-noreason`)
/ Implementer (`deepseek-coder-v2-lite`) loop over the T7610 router, not
authored directly.** An earlier pass through items 2-7 was done by direct
authorship for speed; it was discarded and redone once that was caught,
since testing the local pipeline — not just producing docs — is the
actual point of this project.

## Items 1-7: DONE (22 pages built, matching this scope's file count)

### Category 1 — Direct extraction (core/, item 1 + item 2's main page)

| File | Notes |
|---|---|
| `core/PILOT.m` | First page built; also the §4a toolchain trial's test case |
| `core/SIFTED.m` | One retry (split combined row); one manual fix (misattributed text on wrong field) |
| `core/CLOISTER.m` | Manual fix: added missing parent `opts`/`out` summary rows |
| `core/PYTHIA.m` | Cleanest Implementer output of the batch |
| `core/FILTER.m` | One retry (split combined row); manual fix: restored truncated description |
| `core/PILOTviewpoint.m` | Manual fix: added missing parent rows, restored truncated clauses |
| `core/INIT.m` | Clean on first pass |
| `core/PRELIM.m` | Full manual rewrite after two unreliable Implementer attempts |
| `core/TRACE.m` (item 2) | Implementer invented type annotations (`char --`, `double --`, etc.) not present in the source — stripped directly |

### Category 2 — Docstring expansion required first (utils/, output/, wrappers: item 4-6, 11 files)

Two-stage pipeline: Scaffold extracted grounding facts from each
function's actual code; Planner+Implementer drafted an Inputs/Outputs
docstring addition from those facts (committed to the `.m` source);
Planner+Implementer then extracted the page from the now-complete
docstring, gated the same way as Category 1.

| File | Notes |
|---|---|
| `utils/ISAdefaults.m` | Source bug found+fixed (see below); one duplicate-paragraph cleanup |
| `utils/ISAvalidateOpts.m` | Source bug found+fixed; content was correct both attempts, just unreachable until the source was fixed |
| `utils/ISAsubsetData.m` | Source bug found+fixed; one completeness fix (Description dropped the featIdx clause) |
| `utils/ISAgetClassifierFcn.m` | One formatting nit (missing space in a table cell) |
| `utils/ISAmigrateModel.m` | Source bug found+fixed; Implementer dropped the Inputs/Outputs section twice on this file's long docstring and once "corrected" a documented typo (`bestPerformace`) it should have preserved — completed directly from verified source text after two failed retries, per the framework's bounded-retry-then-Scaffold rule |
| `output/scriptcsv.m` | Draft-stage Planner initially told the Implementer it could omit the Outputs section for a void function — fixed in the shared prompt |
| `output/scriptpng.m` | Clean |
| `output/scriptweb.m` | Clean |
| `output/ISArecallView.m` | Clean |
| `buildIS.m` | Source bug found+fixed; first attempt (on the truncated source) fabricated a nonexistent `model_struct` variable name and invented descriptions |
| `exploreIS.m` | Source bug found+fixed; same fabrication pattern as buildIS.m on the first attempt |

`core/TRACE_legacy.m` (item 2's folded subsection) also went through this
stage: Inputs/Outputs drafted from facts, committed to source, then
synthesised into TRACE.md's "Legacy method" subsection.

### Category 3 — Narrative pages (item 3, item 7)

| Page | Notes |
|---|---|
| `InstanceSpace.m` class overview (item 3) | The hard spot per §4c. First attempt stopped ~40% through a long fact-sheet without hitting its token budget — same "drops trailing content on long documents" failure as the Category 2 long-docstring cases, now confirmed on synthesis too, not just extraction. Fixed by splitting the fact-sheet into two parts (Overview/Constructor/Properties/save/load, then Methods/Usage) and concatenating, with Scaffold merging and de-duplicating the two Implementer outputs |
| Deprecated functions combined page (item 7) | Clean; one completeness fix (dropped the "SIFTED2 was renamed to SIFTED" detail) |

## Real pipeline bugs found and fixed, recorded here and in the framework doc

1. **The Planner must be `gemma-noreason`, not `gemma`.** Reasoning-on
   silently burns the whole completion budget on `reasoning_content` and
   never emits an answer (`finish_reason: "length"`, `content: ""`) —
   every one of item 1's first 7 Planner calls returned empty before this
   was caught.
2. **A genuine source bug, not just a doc-extraction artifact.** The
   docstring-insertion script left a literal blank line between a
   function's original docstring and the newly-drafted Inputs/Outputs
   section in most of the Category 2 files. A blank line ends a MATLAB
   help-comment block — `help`/`doc` on these functions would never have
   shown the new section either, not just this pipeline's extraction.
   Fixed in every affected file by converting the blank line to a bare
   `%` continuation.
3. **The Implementer fabricates content rather than reporting absence.**
   When the source bug above made the real content unreachable
   (`buildIS`/`exploreIS`), the Implementer invented a plausible-looking
   `model_struct` variable name and paraphrased descriptions instead of
   flagging that nothing was there. Confirmed again independently on
   `TRACE.m`, where it invented MATLAB type annotations (`char --`,
   `double --`) matching the *style* of other files' opts tables but not
   present in TRACE.m's own source.
4. **The Implementer drops trailing sections on sufficiently long
   documents, even when the source is intact and well within its token
   budget** (`finish_reason` was `"stop"`, not `"length"`, in every case
   checked). Confirmed three times: `ISAvalidateOpts.m`/`ISAmigrateModel.m`
   twice each (docstring extraction) and `InstanceSpace.md` once
   (narrative synthesis). An explicit "don't stop early" instruction
   fixed nothing; splitting the input into two shorter passes did.

## Items 8, 9, 12, 13: DONE

| Page | Notes |
|---|---|
| Landing (item 8) | Clean on first pass -- short, well-scoped facts avoided the long-document issue entirely |
| GettingStarted (item 9) | Clean; code block verified to match `example.m` exactly |
| MetadataFormat (item 12) | The whole schema was written out twice (once cleanly, once as a near-duplicate restatement) -- trimmed the duplicate half |
| MigratingLegacyModel (item 13) | Clean; fixed the cross-link to `ISAmigrateModel.html` to match the site's link convention |

## Item 10: DONE

Interactive Walkthrough (`liveDemoIS.m`), split into two parts up front
given the demonstrated long-document risk -- both came back clean and
accurate on the first pass. Two small fixes: both parts independently
used LaTeX-style math notation (`\(...\)`, `$Z$`) that won't render on
this plain-Markdown site, converted to inline code (`` `Z = X*A'` ``);
and Part 1's own closing sentence had to be dropped after merging, since
it read as a premature conclusion sitting in the middle of the combined
page.

## Item 11: DONE

Options Reference, split into two parts (general/perf/prelim/auto/bound/
norm/selvars, then sifted/pilot/cloister/pythia/trace/outputs), every
field/default checked one by one against `utils/ISAdefaults.m`. Two real
findings:
1. **Part 1 invented "feature selection" framing for `opts.selvars`'s
   under-specified fields** (`smallscaleflag`, `fileidxflag`,
   `densityflag`, etc.) — plausible-sounding, but wrong. Verified against
   `InstanceSpace.m` directly (`cvpartition`/file-index/FILTER-density
   logic around line 580): this group is **instance subsetting**, not
   feature selection (that's SIFTED's job, a different stage entirely).
   Corrected all six descriptions from the real source.
2. **Part 2 truncated again on a long table** — missing `opts.outputs`'s
   `fig`/`web` rows and a few empty description cells at the very end.
   Completed directly from already-verified facts rather than a third
   retry of the same prompt.
3. **A misleading default value**: `opts.trace.contra`'s table showed a
   bare `true`, but the actual default (given `method='trace3'` is the
   real default) is `false` — `true` only applies when
   `method='legacy'`. Reworded to `false (true when method='legacy')`.

## Item 14: DONE

What's New, deliberately condensed per the issue's own wording ("can
largely point at RELEASE_NOTES.md, or be its own condensed page") —
headline bullets only for v0.9.1/v0.9.0, pointing to RELEASE_NOTES.md
for the full list. Two fixes: the recurring `ISAmigrateModel` naming
typo (wrong capitalization, seen before in item 3's draft) appeared
twice; and the closing reference link was a **fabricated placeholder
GitHub URL** (`https://github.com/your-repo/RELEASE_NOTES.md`) that
doesn't point anywhere real — replaced with a plain-text file reference.

**All 29 content pages are now complete, exactly matching issue #53's
own stated inventory (22 reference + 7 conceptual).**

## Items 15-16: DONE — the entire 16-item registry is complete

| File | Notes |
|---|---|
| `info.xml` (item 15) | Schema independently verified against MathWorks documentation (not just trusted from the framework doc's own summary, which omitted the required `<icon>` element) — `<productinfo>` root with `matlabrelease`/`name`/`type`/`icon`/`help_location`/`help_contents_icon` in that exact order. Scaffold caught and fixed one real error the facts themselves introduced: `help_location` was given as bare `html`, but `info.xml` sits at the repo root while the real HTML lives at `doc/html/` — corrected to `doc/html` |
| `doc/html/helptoc.xml` (item 15) | All 29 pages present, none missing, none invented (verified programmatically by parsing the XML and diffing target sets). **Real structural deviation**: the facts asked for 9 grouped sections mirroring `Contents.m`'s categories (an explicit design requirement, not a nice-to-have) — the Implementer flattened almost everything to one top-level list, keeping only one section properly nested. Rebuilt directly with the correct nesting. A second, self-caught error: the first rebuild nested all 9 sections *inside* the Landing page's own `tocitem` instead of as its siblings, contradicting the facts' own "first tocitem is the landing page" framing (and the MathWorks convention it's based on) — fixed before landing |
| `.github/workflows/docs-pages.yml` (item 16) | Clean and accurate on the first pass — correctly adapted from `pyInstanceSpace`'s real workflow (no Python/pdoc/poetry steps, since this repo's `doc/html/` is pre-built static HTML, not CI-generated; correct branch `master`; correct `path: doc/html`). Verified as syntactically valid YAML (the `on:`-parses-as-boolean quirk in PyYAML's `safe_load` is a well-known non-issue, not specific to this file — GitHub's own parser handles it correctly) |

## Summary: the entire 16-item task registry is complete

29 content pages (22 reference + 7 conceptual, exactly matching issue
#53's own stated inventory) + `info.xml` + `helptoc.xml` +
`docs-pages.yml`. Every item went through the real Planner/Implementer
loop over the T7610 router, with Scaffold verification catching a real,
distinct problem in the majority of items — this was as much a stress
test of the local pipeline as it was a documentation build, per the
project's actual objective.

**Full accumulated list of real pipeline bugs found and fixed across the
whole project:**
1. Planner (`gemma`) reasoning mode silently returning empty output —
   fixed by switching to `gemma-noreason`.
2. A genuine MATLAB source bug (not just a doc-extraction artifact): a
   blank line breaking the comment block, found in most Category 2
   files.
3. The Implementer fabricating plausible content (a nonexistent
   `model_struct` variable, invented MATLAB type annotations) when given
   truncated or under-specified input, rather than reporting absence.
4. The Implementer dropping trailing content on long documents even well
   within its token budget (`finish_reason: "stop"`) — confirmed on both
   extraction and narrative synthesis; fixed by splitting input, not by
   asking it to try harder.
5. An invented "feature selection" framing for `opts.selvars` fields
   that are actually about instance subsetting — caught by checking the
   real source (`InstanceSpace.m`), not by trusting a plausible
   description.
6. A misleading conditional default value (`opts.trace.contra`) stated
   as a bare `true` when the real default (under `method='trace3'`) is
   `false`.
7. A fabricated placeholder URL (a fake GitHub link) in the What's New
   page's closing reference.
8. A structural spec deviation in `helptoc.xml` — flattening a required
   nested category structure — caught by checking against the facts'
   own explicit design, not just validating well-formedness.

## Fix: cross-page navigation (2026-09-10)

Confirmed gap: 20 of 29 pages had zero links to any other page; the
other 9 had only incidental prose links; `generate.py`'s `TEMPLATE` had
no nav markup at all. Fixed per
`llm-server/docs/instancespace_docs_nav_framework.md` — a shared
sidebar, sourced from `helptoc.xml`'s own 29-target taxonomy (not
re-derived), injected into every page via a new `NAV_DATA` constant and
`render_nav()`. Executed through the real local pipeline
(`gemma-noreason` Planner, `qwen3-coder-next` Implementer, both used as
already-loaded, no reload). Passed the mechanical gate on the first
attempt — no bounded retries needed, a first in this project's history.
Sidebar CSS (flex layout, `.current-page` marker) was added directly by
Scaffold afterward — the Implementer's HTML/data was correct but had no
corresponding stylesheet rule, so it would have rendered as two stacked
blocks rather than an actual sidebar; caught by an actual browser
screenshot, not by the mechanical gate, which only checked link
correctness, not layout. Verified end-to-end: screenshot + a real click
from PILOT → SIFTED confirmed the sidebar renders and every link
navigates correctly.

## Fix: content defects surfaced by actually browsing the site (2026-09-10)

The user reported 6 defects after browsing the nav-fixed site by hand —
none introduced by the nav fix, all pre-existing in the original 16-item
build despite it being recorded as complete. All verified directly
before fixing (per this framework's own "audit before you tier"
principle) rather than trusted at face value:

1. **10 pages had a literal `TODO: cross-link related functions` stub**
   under `## See Also`, never filled in: `InstanceSpace`, `PRELIM`,
   `ISAmigrateModel`, `PYTHIA`, `INIT`, `SIFTED`, `ISAsubsetData`,
   `TRACE`, `OptionsReference`, `InteractiveWalkthrough`.
2. **Landing's title was literally "Landing"** and carried a
   hand-authored, non-clickable "Documentation Navigation" bullet list
   duplicating the sidebar — removed; title corrected to "InstanceSpace"
   per `helptoc.xml`'s own approved label for that page.
3. **Copyright footer** existed in only 2 of 29 pages (inconsistent
   because it was per-source-file text, not template-driven) — moved
   into `generate.py`'s `TEMPLATE` so all 29 get it uniformly; the 2
   duplicate inline copies removed.
4. **`CLOISTER.md`/`FILTER.md` had no See Also, References, or
   Copyright at all** — ended right after Output Arguments, unlike
   every other function page.
5. **`PILOTviewpoint.md`'s math notation
   (`||Y_group - (C*A*Z')'||_F^2 + ...`) was plain prose text**, not
   wrapped in a code span the way the identical situation in `PILOT.md`
   already was — fixed to match that existing convention (no MathJax
   added; not warranted for this).
6. **References existed on only 3 of 29 pages** — every function-backed
   `.m` source file already carries its own verbatim `Reference:`
   comment block; extracted and formatted to match `PILOT.md`'s
   established citation style for the 18 pages missing it, via a
   deterministic script (zero LLM involvement — pure extraction, no
   content invented).

**A 7th, previously unreported defect found while fixing #6**:
`scriptcsv.md` and `ISArecallView.md` had **wrong** content — both had
`PILOT.md`'s own References/See Also copy-pasted onto them verbatim
(SIFTED/PILOTviewpoint links on a CSV-writer function has nothing to do
with what it does). Caught only because the extraction script's
"already has References" skip made this an active investigation instead
of silently trusting existing content. `scriptcsv.md`'s References
corrected to its own single real citation; both files' See Also
corrected to their real relationships (confirmed via `InstanceSpace.m`
call-graph grep, not guessed).

See Also content for all 14 affected pages (10 TODO + `CLOISTER` +
`FILTER` + the 2 contaminated files) was drafted by the Implementer
(`qwen3-coder-next`) from real, grep-verified relationships (pipeline
call order and specific call-site evidence in `InstanceSpace.m`, e.g.
"SIFTED calls PILOT internally," "PYTHIA's Yhat feeds TRACE") — the
Implementer was given the facts and told explicitly not to invent
beyond them; all 14 outputs matched the given facts exactly with zero
extra links, a first-shot pass with no corrections needed.

**Known, deliberately out-of-scope observation**: `ISAvalidateOpts`,
`ISAgetClassifierFcn`, `scriptpng`, `scriptweb`, `Deprecated`,
`MetadataFormat`, and `WhatsNew` also have no See Also section — but
none of these were ever a `TODO` stub (they simply never had one), and
none were part of the user's reported list. Not touched, to avoid scope
creep beyond what was actually reported.

## Current repo state

Branch `docs/site-skeleton` in `andremun/InstanceSpace`, all changes
uncommitted. Modified: 12 `.m` source files (11 Category-2
Inputs/Outputs additions + `TRACE_legacy.m`), `doc/generate.py` and
`doc/style.css` (nav fix above). Added: `doc/` (29 `.md` + 29 `.html` +
`helptoc.xml` + `style.css` + `generate.py`), `info.xml`,
`.github/workflows/docs-pages.yml`, this status file. Nothing has been
committed or pushed — ready for review.
