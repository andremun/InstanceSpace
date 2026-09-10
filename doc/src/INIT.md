# INIT
Load and filter instance-space metadata (build or explore time).
## Syntax
```
Training mode (2 args): [data, extra] = INIT(rootdir, opts)
Evaluation mode (3 args): [data, extra] = INIT(rootdir, opts, trainedModel)
```
## Description
Reads rootdir/metadata.csv, applies opts.selvars.feats/.algos column filtering, and drops instances/features with too many missing values (opts.prelim.nanThreshold).

Reads rootdir/metadata_test.csv, applies the same opts.selvars filtering, validates the feature set/order against trainedModel (a previously-built InstanceSpace object's obj.model), and reconciles algorithm columns against trainedModel.data.algolabels: known algorithms line up by name at their trained column position, unseen ones are appended as new columns (NaN for training-only algorithms).

Before this function existed, this logic was two independent, drifted implementations -- the version embedded in InstanceSpace.runPrelim (build time) and a separate inline reimplementation in InstanceSpace.evaluateTestSet (explore time). INIT is a pure extraction of both into one shared function, dispatched by nargin like PYTHIA/TRACE (#38): each mode's behaviour is preserved exactly as it was, just no longer duplicated in two places that could silently drift apart from each other.
## Input Arguments
| Argument | Description |
|---|---|
| rootdir | directory containing metadata.csv (training mode) or metadata_test.csv (evaluation mode), trailing slash |
| opts | the full opts struct (obj.opts in training mode, trainedModel.opts in evaluation mode) |
| trainedModel | (optional) a previously-built InstanceSpace object's obj.model; presence selects evaluation mode |
## Output Arguments
| Argument | Description |
|---|---|
| data | struct: instlabels, S (if a 'source' column is present), X, Y, Xraw, Yraw, algolabels. Training mode also: featlabels (evaluation mode's equivalent is finalised later by the caller, after featsel.idx subsetting -- see extra.featlabelsAll) |
| extra | training mode: struct() (unused) evaluation mode: struct with featlabelsAll (feature labels, 'feature_'-prefixed, after opts.selvars.feats filtering but before featsel.idx subsetting), modelalgos, newalgos, nalgos, ninst -- all needed by the caller's downstream PRELIM/normalisation/featsel step, which must treat the first modelalgos columns of data.Y differently from newly appended ones |
## References

- Smith-Miles, K. & Munoz, M.A. (2023). Instance Space Analysis for Algorithm Testing. ACM Computing Surveys, 55(12), Article 255. <https://doi.org/10.1145/3572895>

## See Also

[PRELIM](PRELIM.html) | [InstanceSpace](InstanceSpace.html)
