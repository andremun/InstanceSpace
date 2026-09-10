# ISAvalidateOpts

Validate user-supplied opts fields before defaults are filled.

## Syntax

```
opts = ISAvalidateOpts(opts)
```

## Description

checks the type/range of every RECOGNISED opts field the caller actually supplied (the fixed set of fields this function knows about; see the body) and errors clearly (ISA:ISAvalidateOpts:*) on the first invalid one, instead of letting an out-of-range value surface many stages later as a confusing crash deep inside PRELIM/PILOT/PYTHIA/etc. An unrecognised field name (e.g. a typo like opts.piyhia.classifier) is NOT flagged -- it passes through silently, exactly like an unset one, since this function has no way to distinguish "not a real option" from "a future option it doesn't know about yet". Deliberately validates only fields that ARE present: this runs before ISAdefaults, so most fields are still absent at this point and are not this function's concern -- ISAdefaults supplies known-valid defaults for anything missing. "Present" is tracked explicitly (getf() returns a presence flag alongside the value), not inferred from isempty(v): a field explicitly supplied as opts.general.parallel = [] IS present and must be rejected as invalid, not silently skipped as if absent -- ISAdefaults only checks isfield(), so it would never replace that [] with the proper default, and later code expecting a logical scalar would fail far from the actual mistake. opts is returned unmodified; this function only ever errors or passes through, it never rewrites values (renaming/migrating legacy field names is ISAmigrateModel's job, not this one's).

## Input Arguments

| Argument | Description |
|---|---|
| `opts` | struct, user-supplied, possibly partial -- only the fields actually present are checked |

## Output Arguments

| Argument | Description |
|---|---|
| `opts` | the same struct, unmodified. This function only ever errors (on the first invalid recognised field) or passes through -- it never rewrites values |

## References

- Smith-Miles, K. & Munoz, M.A. (2023). Instance Space Analysis for Algorithm Testing. ACM Computing Surveys, 55(12), Article 255. <https://doi.org/10.1145/3572895>
