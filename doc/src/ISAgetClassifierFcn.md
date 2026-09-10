# ISAgetClassifierFcn

Resolve a classifier name from the registry.

```
[fitFcn, p1label, p2label] = ISAgetClassifierFcn(name)
```

## Registry:

| name      | MATLAB fn     | param 1             | param 2          |
|-----------|---------------|---------------------|------------------|
| 'knn'     | fitcknn       | NumNeighbors [1,25] | Distance (cat.)  |
| 'svm'     | fitcsvm       | BoxConstraint log2 | KernelScale log2 |
| 'tree'    | fitctree      | MinLeafSize [1,100] | N/A              |
| 'nb'      | fitcnb        | Bandwidth  log10   | N/A              |
| 'linear'  | fitclinear    | Lambda     log10   | N/A              |
| 'ensemble'| fitcensemble  | NumLearningCycles [10,200] | MinLeafSize [1,20]|

fitcecoc is excluded: PYTHIA trains one binary classifier per algorithm; multi-class ECOC machinery is never required.

## Input Arguments

| Argument | Description |
|----------|-------------|
| `name`   | string representing the classifier registry entry (options: 'knn', 'svm', 'tree', 'nb', 'linear', 'ensemble'). |

## Output Arguments

| Field | Description |
|-------|-------------|
| `fitFcn` | MATLAB fitc* function handle. |
| `p1label` | human-readable label for the first Sobol-tuned hyperparameter. |
| `p2label` | human-readable label for the second Sobol-tuned hyperparameter (may be 'N/A' for classifiers with only one tunable parameter). |

## References

- Smith-Miles, K. & Munoz, M.A. (2023). Instance Space Analysis for Algorithm Testing. ACM Computing Surveys, 55(12), Article 255. <https://doi.org/10.1145/3572895>
