# Metadata File Format

Prepare metadata.csv for your own problem domain

An instance space is built from one CSV file, `metadata.csv`, with one row per problem instance. Instances to evaluate later with `explore` go in `metadata_test.csv`, in the same format.

## Columns

Columns are recognised by their header. The order of the column groups does not matter, and headers are case insensitive.

| Header | Required | Contents |
|---|---|---|
| `Instances` | yes | Instance name. Text or numbers. |
| `Source` | no | Where the instance comes from, for example the benchmark suite. Drawn by `obj.plot('sources')`. |
| `feature_<name>` | at least two | Numeric feature value. |
| `algo_<name>` | at least one | Numeric performance of the algorithm on the instance. |

Other columns are ignored. Names after the prefix must be valid MATLAB identifiers: use underscores instead of spaces, for example `feature_edge_density` or `algo_simulated_annealing`.

## Example

The first rows of the reference data (`test/data/metadata.csv`), shortened to three features and three algorithms:

```
Instances,feature_Max_Normalized_Entropy_attributes,feature_ErrorRate_Decision_Node,feature_Training_Error_Linear_Classifier_L2,algo_NB,algo_KNN,algo_RandF
abalone,0.332548382,0.386458668,0.260234614,0.280095763,0.256969153,0.234432252
abalone_ori,0.332548382,0.210409617,0.035714286,0.164181841,0.167112964,0.895471859
```

## Missing Values

Leave a cell empty or write `NaN` for a missing value. A feature with at least `opts.prelim.nanThreshold` missing values is removed. A missing performance value never counts as good.

Avoid `NA`, spreadsheet error codes such as `#DIV/0!`, and empty rows: they make MATLAB read the whole column as text, which then fails later in the pipeline.

## Test Instances

`metadata_test.csv` must contain the same `feature_` columns as `metadata.csv`, in the same order. It may contain:

- a subset of the trained algorithms — the missing ones are not scored;
- algorithms that were not trained — they are appended to the results, without a classifier.

## Choosing Features

The features determine what the instance space can reveal. Good features:

- are cheap to compute compared with running the algorithms;
- capture properties that plausibly affect how hard an instance is for each algorithm;
- vary across the instances, with few repeated values.

Start with more candidate features than you need; `SIFTED` removes the ones that do not explain performance. See Smith-Miles & Muñoz (2023) for guidance on feature design.

## References

- Smith-Miles, K. & Muñoz, M.A. (2023). Instance Space Analysis for Algorithm Testing. *ACM Computing Surveys*, 55(12), Article 255. <https://doi.org/10.1145/3572895>

## See Also

`INIT` | `InstanceSpace` | [Options Reference](OptionsReference.html)
