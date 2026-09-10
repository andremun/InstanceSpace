# exploreIS

Thin backward-compatibility wrapper around InstanceSpace.

Requires model.mat to already exist in rootdir (written by buildIS). Preserves the original exploreIS(rootdir) calling convention for callers -- notably the MATILDA web platform -- that invoke this entry point directly. New code should use InstanceSpace directly:

```
obj = InstanceSpace.load(rootdir);
obj = obj.explore(rootdir);
out  = obj.getResults(1);
```

## Input Arguments

| Argument | Description |
|---|---|
| `rootdir` | string/path to the directory containing model.mat and metadata_test.csv |

## Output Arguments

| Field | Description |
|---|---|
| `out` | the most recent test-results entry (specifically obj.testResults{end}) |

## References

- Smith-Miles, K. & Munoz, M.A. (2023). Instance Space Analysis for Algorithm Testing. ACM Computing Surveys, 55(12), Article 255. <https://doi.org/10.1145/3572895>

## See Also

[InstanceSpace](InstanceSpace.html)
