# ISAdefaults

Fill in missing opts fields with default values.

## Syntax

```
opts = ISAdefaults(opts)
```

## Description

ensures every pipeline function receives a complete options struct, eliminating scattered isfield chains across buildIS, PILOT, TRACE, and CLOISTER. Call once at the buildIS entry point after loading options.json.

## Input Arguments

| Argument | Description |
|---|---|
| `opts` | struct, potentially partial (e.g., freshly parsed from options.json with some fields absent or left as []) |

## Output Arguments

| Field | Description |
|---|---|
| `opts` | the same struct, with every pipeline-stage field guaranteed present. See the Options Reference page for the full field list and defaults. |

## References

- Smith-Miles, K. & Munoz, M.A. (2023). Instance Space Analysis for Algorithm Testing. ACM Computing Surveys, 55(12), Article 255. <https://doi.org/10.1145/3572895>

## See Also

[buildIS](buildIS.html) | [PILOT](PILOT.html) | [TRACE](TRACE.html) | [CLOISTER](CLOISTER.html)
