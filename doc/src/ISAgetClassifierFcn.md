# ISAgetClassifierFcn

Look up a classifier in the PYTHIA registry

## Syntax

```
[fitFcn,p1label,p2label] = ISAgetClassifierFcn(name)
```

## Description

`[fitFcn,p1label,p2label] = ISAgetClassifierFcn(name)` returns the Statistics and Machine Learning Toolbox fitting function for the classifier `name`, and the names of the hyperparameters that `PYTHIA` tunes for it. `name` is the value of `opts.pythia.classifier`.

| `name` | Fitting function | First hyperparameter (search range) | Second hyperparameter |
|---|---|---|---|
| `'knn'` | `fitcknn` | `NumNeighbors` [1, 25] | `Distance` (categorical) |
| `'svm'` | `fitcsvm` | `BoxConstraint` (log2 scale) | `KernelScale` (log2 scale) |
| `'tree'` | `fitctree` | `MinLeafSize` [1, 100] | — |
| `'nb'` | `fitcnb` | `Bandwidth` (log10 scale) | — |
| `'linear'` | `fitclinear` | `Lambda` (log10 scale) | — |
| `'ensemble'` | `fitcensemble` | `NumLearningCycles` [10, 200] | `MinLeafSize` [1, 20] |

Each classifier is binary: PYTHIA trains one per algorithm, so multiclass learners such as `fitcecoc` are not needed.

## Examples

### Look up the SVM entry

```matlab
[fitFcn, p1, p2] = ISAgetClassifierFcn('svm')
```

```
fitFcn = @fitcsvm
p1 = 'BoxConstraint'
p2 = 'KernelScale'
```

### Supply fixed hyperparameters for a single-parameter classifier

With `opts.pythia.tuning = 'none'`, `opts.pythia.params` needs one column per hyperparameter.

```matlab
[~, ~, p2] = ISAgetClassifierFcn('tree');
nparams = 1 + ~strcmp(p2, 'N/A');     % 1 for a tree
opts.pythia.classifier = 'tree';
opts.pythia.tuning = 'none';
opts.pythia.params = repmat(10, 10, nparams);   % MinLeafSize 10 for 10 algorithms
```

## Input Arguments

### `name` — Classifier name

*`'knn'` | `'svm'` | `'tree'` | `'nb'` | `'linear'` | `'ensemble'`*

Case insensitive. Any other value raises `ISA:ISAgetClassifierFcn:unknownClassifier`.

## Output Arguments

### `fitFcn` — Fitting function

*function handle*

### `p1label` — First hyperparameter name

*character vector*

### `p2label` — Second hyperparameter name

*character vector*

`'N/A'` for classifiers with one hyperparameter.

## Version History

### v0.9.0 — Introduced

## See Also

`PYTHIA` | [Options Reference](OptionsReference.html#opts-pythia)
