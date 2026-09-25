# Instance Space Analysis: A toolkit for the assessment of algorithmic power

[![View InstanceSpace on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://au.mathworks.com/matlabcentral/fileexchange/75170-instancespace)
[![DOI](https://zenodo.org/badge/144672744.svg)](https://doi.org/10.5281/zenodo.4484107)
[![Citations](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fapi.semanticscholar.org%2Fgraph%2Fv1%2Fpaper%2FDOI%3A10.1145%2F3572895%3Ffields%3DcitationCount&query=%24.citationCount&label=citations&color=blue)](https://doi.org/10.1145/3572895)
[![Downloads](https://img.shields.io/github/downloads/andremun/InstanceSpace/total.svg)](https://github.com/andremun/InstanceSpace/releases)
[![Tests](https://github.com/andremun/InstanceSpace/actions/workflows/tests.yml/badge.svg)](https://github.com/andremun/InstanceSpace/actions/workflows/tests.yml)
[![codecov](https://codecov.io/gh/andremun/InstanceSpace/branch/master/graph/badge.svg)](https://codecov.io/gh/andremun/InstanceSpace)
[![Documentation](https://img.shields.io/badge/docs-reference-blue)](https://andremun.github.io/InstanceSpace/)

Instance Space Analysis is a methodology for assessing the strengths and weaknesses of an algorithm and objectively comparing its algorithmic power, without bias introduced by a restricted choice of test instances. At its core is the modelling of the relationship between an instance's structural properties and the performance of a group of algorithms. Instance Space Analysis allows the construction of **footprints** for each algorithm, defined as regions in the instance space where we statistically infer good performance. Other insights that can be gathered from Instance Space Analysis include:

-	Objective metrics of each algorithm’s footprint across the instance space as a measure of algorithmic power;
-	Explanation through visualisation of how instance features correlate with algorithm performance in various regions of the instance space;
-	Visualisation of the distribution and diversity of existing benchmark and real-world instances;
-	Assessment of the adequacy of the features used to characterise an instance;
-	Partitioning of the instance space into recommended regions for automated algorithm selection;
-	Distinguishing areas of the instance space where it may be useful to generate additional instances to gain further insights.

The unique advantage of visualising algorithm performance in the instance space, rather than as a small set of summary statistics averaged across a selected collection of instances, is the nuanced analysis it enables: explaining strengths and weaknesses and examining interesting variations in performance that may be hidden by tables of summary statistics.

This repository provides a set of MATLAB tools for performing a complete Instance Space Analysis as part of an automated pipeline. It is also the computational engine that powers the Melbourne Algorithm Test Instance Library with Data Analytics ([MATILDA](http://matilda.unimelb.edu.au/matilda/)) web tools for online analysis. For further information on the Instance Space Analysis methodology, see [here](http://matilda.unimelb.edu.au/matilda/our-methodology).

If you follow the Instance Space Analysis methodology, please cite as follows:

> K. Smith-Miles and M.A. Muñoz. *Instance Space Analysis for Algorithm Testing: Methodology and Software Tools*. ACM Comput. Surv. 55(12:255),1-31 [DOI:10.1145/3572895](https://doi.org/10.1145/3572895), 2023.

If you use the 3D extension of the methodology (`opts.pilot.dims = 3`, `PILOTviewpoint`, TRACE3's native 3D footprints), please additionally cite:

> C. Simpson, M.A. Muñoz, S. Kandanaarachchi and R.J.G.B. Campello. *ISA3: A 3-dimensional expansion of Instance Space Analysis*. Machine Learning, 114, 240 [DOI:10.1007/s10994-025-06871-5](https://doi.org/10.1007/s10994-025-06871-5), 2025.

Also, if you specifically use this code, please cite as follows:

> M.A. Muñoz and K. Smith-Miles. *Instance Space Analysis: A toolkit for the assessment of algorithmic power*. andremun/InstanceSpace on GitHub. Zenodo, [DOI:10.5281/zenodo.4484107](https://doi.org/10.5281/zenodo.4484107), 2020.

Or if you specifically use [MATILDA](http://matilda.unimelb.edu.au/matilda/), please cite as follows:

> K. Smith-Miles, M.A. Muñoz and Neelofar. *Melbourne Algorithm Test Instance Library with Data Analytics (MATILDA)*. Available at (https://matilda.unimelb.edu.au). 2020.

**DISCLAIMER: This repository contains research code. On occasion, new features will be added, or changes will be made that may result in crashes. Although we have made every effort to reduce bugs, this code comes with NO GUARANTEES. If you find any issues, let us know ASAP via the contact methods listed at the end of this document.**

## Documentation

The reference documentation is at **<https://andremun.github.io/InstanceSpace/>**: a page for every function and for the `InstanceSpace` class, with examples, a getting-started guide, a stage-by-stage walkthrough, the metadata file format, and the full options reference. The same pages open inside MATLAB: with the repository root on the path, open the Help browser and choose **Supplemental Software → Instance Space Analysis Toolbox**.

## Installation Instructions

The main requirement for the software to run is to have MATLAB R2025a or later, with the [Global Optimization](https://au.mathworks.com/help/gads/index.html), [Parallel Computing](https://www.mathworks.com/products/parallel-computing.html), [Optimization](https://au.mathworks.com/products/optimization.html), [Statistics and Machine Learning](https://au.mathworks.com/help/stats/index.html), and [Financial](https://au.mathworks.com/products/finance.html) toolboxes installed. The LIBSVM MEX-files (`svmpredict`/`svmtrain`), used for legacy models, are **not bundled with this repository** (no build source is available for them). They're only relevant for evaluating a pre-v0.9.0 model whose classifiers `ISAmigrateModel` couldn't retrain (missing original training data); `ISAmigrateModel` prefers retraining from scratch with the current registry whenever the training data is available, which needs no LIBSVM dependency. If you do hit that path, obtain LIBSVM yourself from [the official project](https://www.csie.ntu.edu.tw/~cjlin/libsvm/) and add its MEX-files to the MATLAB path — `PYTHIA`'s eval mode raises a clear `ISA:PYTHIA:noLibsvm` error naming the missing dependency if you don't.

## Repository layout

```
InstanceSpace.m, buildIS.m, exploreIS.m   entry points (see below)
example.m, test_integration.m             getting-started / regression suite
liveDemoIS.m                              interactive, stage-by-stage walkthrough (open in
                                          MATLAB's Live Editor)
startup.m                                 adds the folders below to the MATLAB path
Contents.m                                MATLAB Central File Exchange version/date metadata
CITATION.cff                              machine-readable citation metadata
core/                                     PRELIM, SIFTED, PILOT, PILOTviewpoint,
                                          CLOISTER, PYTHIA, TRACE, TRACE_legacy, FILTER
output/                                   scriptcsv, scriptpng, scriptweb, scriptfcn,
                                          ISArecallView
utils/                                    ISAdefaults, ISAvalidateOpts, ISAgetClassifierFcn,
                                          ISAmigrateModel, ISAsubsetData
deprecated/                               PYTHIA2, PYTHIAtest, SIFTED2 (warn-and-forward
                                          shims kept for backward compatibility)
doc/                                      reference documentation: Markdown sources (src/),
                                          generated pages (html/), and generate.py
info.xml                                  registers doc/html with the MATLAB Help browser
```

`InstanceSpace.m`/`buildIS.m`/`exploreIS.m` add `core/`, `output/`, `utils/`, and `deprecated/` to the MATLAB path automatically the first time any of them is used in a session — `example.m`, `test_integration.m`, and any script that starts with `buildIS`/`exploreIS`/`InstanceSpace` need no extra setup. If you want to call a function from one of those folders directly (`PILOT`, `SIFTED`, ...) without going through one of those three first, run `startup.m` (e.g. `run('startup.m')`, or just `startup` with the repo root as your current folder) at the start of your session.

## Working with the code

Start with `example.m`: it runs the full pipeline (`buildIS` + `exploreIS`) once, on the bundled reference dataset, with sensible defaults and just a handful of commonly adjusted settings (classifier, tuning strategy, projection dimensionality, feature selection on/off) exposed as plain variables near the top. Outputs — images (`.png`), tables (`.csv`), and raw intermediate data (`.mat`) — land in `test/data/example/`. To analyse your own data, point it at a folder containing your `metadata.csv` instead (see "The metadata file" below), and revisit the performance-metric settings, which are tuned for the bundled dataset's error-rate semantics.

`test_integration.m` is a thin `matlab.unittest` runner over the exhaustive option-coverage regression suite in `tests/*.m` (every classifier, tuning strategy, 2D/3D, PLS, viewpoint groups, staged `build()`/`explore()`/save-load round-trips, and the full `ISAmigrateModel` legacy-migration table), each case in its own subdirectory under `test/data/` (e.g. `test/data/classifier_svm/`) so no run overwrites another's outputs, with a code-coverage report (`coverage.xml`) produced alongside. It's a good reference for how a given option is meant to be used, but not the place to start.

**`options.json` is a generated artifact, not a source file**, for both scripts above. Each run writes its own `options.json` from the `opts` struct built in MATLAB. Hand-editing an `options.json` file has no lasting effect — the next run silently overwrites it. To change what gets run, edit the MATLAB script instead (`example.m` directly, or for the `tests/*.m` suite behind `test_integration.m`, `testDefaultOpts()` for shared settings and a specific `TestParameter` case's `.override` for that case only, e.g. in `PipelineOptionsTest.m`'s `pipelineOptionCases()`).

### The InstanceSpace class

`buildIS`/`exploreIS` are thin backward-compatibility wrappers (kept for callers like the MATILDA web platform that invoke them directly) around the `InstanceSpace` class, which is the recommended interface for new code:

```matlab
obj = InstanceSpace(rootdir);              % reads options.json if present, else defaults
obj = obj.build();                          % run the full pipeline
obj = obj.explore(testRootDir);              % evaluate a trained model on new data
results = obj.getResults();                  % training results (== obj.model)
obj.save();                                  % write rootdir/model.mat
obj = InstanceSpace.load(rootdir);           % read it back
```

Options can be changed between individual pipeline stages, and only the stages that need to re-run do:

```matlab
obj = InstanceSpace(rootdir);
obj = obj.build('stages', {'prelim', 'sifted', 'pilot'});
obj.opts.pilot.alpha = 2.0;                  % adjust after inspecting the projection
obj = obj.build('stages', {'pilot'});         % re-runs PILOT only; sifted output is reused
obj = obj.build('stages', {'cloister', 'pythia', 'trace'});
```

Both `build()` and `explore()` also accept an optional `'onStage'` callback, invoked once after each stage completes with that stage's name and the model/result at that point — useful for inspecting intermediate results (e.g. PILOT's projection) without splitting a run into several staged calls:

```matlab
inspect = @(stageName, model) fprintf('%s done\n', stageName);
obj = obj.build('onStage', inspect);
obj = obj.explore(testRootDir, 'onStage', inspect);
```

See the class's own help text (`help InstanceSpace`) for the full method list, including `plot()` and `getResults(idx)` for accessing a specific `explore()` call's results.

## The metadata file

The `metadata.csv` file should contain a table where each row corresponds to a problem instance, and each column must strictly follow the naming convention mentioned below:

-	**instances** instance identifier - We expect the instance identifier to be of type "String". This column is mandatory.
-	**source** instance source - This column is optional
-	**feature_name** The keyword "feature_" concatenated with feature name. For instance, if the feature name is "density", the header name should be mentioned as "feature_density". If the name consists of more than one word, each word should be separated by "_" (spaces are not allowed). There must be more than two features for the software to work. We expect the features to be of the type "Double".
-	**algo_name** The keyword "algo_" concatenated with algorithm name. For instance, if the algorithm name is "Greedy", the column header should be "algo_greedy". If the name consists of more than one word, each word should be separated by "_" (spaces are not allowed). You can add the performance of more than one algorithm in the same `.csv`. We expect the algorithm performance to be of the type "Double".

Moreover, empty cells, NaN or null values are allowed but **not recommended**. We expect you to handle missing values in your data before processing. You may use [this file](https://matilda.unimelb.edu.au/matilda/matildadata/graph_coloring_problem/metadata/metadata.csv) as reference.

**Common data-preparation mistake**: using `NA` instead of `NaN`, or leaving Excel error codes (`#REF!`, `#NULL!`, `#DIV/0!`) or empty rows in the sheet. Any of these causes `readtable` to infer a column as text (`string`/`cell`) instead of numeric (`double`), which will crash the pipeline downstream rather than failing with a clear error at load time.

## Options

Every setting is a field of the `opts` structure, given to `InstanceSpace` directly or as `options.json` in the data folder. Set only the fields you want to change; every other field takes its default. The [Options Reference](https://andremun.github.io/InstanceSpace/OptionsReference.html) lists every field with its default and meaning. The settings most often changed are:

| Field | Default | Meaning |
|---|---|---|
| `opts.perf.MaxPerf` | `false` | `true` if larger performance values are better; `false` for a cost such as error or run time. |
| `opts.perf.AbsPerf` | `false` | `true`: an algorithm is good when its performance is better than `epsilon`. `false`: when it is within a fraction `epsilon` of the best algorithm on the instance. |
| `opts.perf.epsilon` | `0.05` | Good-performance threshold. |
| `opts.perf.betaThreshold` | `0.55` | An instance is easy when more than this fraction of the algorithms are good on it. |
| `opts.general.seed` | `42` | Random seed for every stochastic stage. |
| `opts.general.parallel` | `false` | Use a parallel pool. |
| `opts.selvars.feats`, `opts.selvars.algos` | all | Cell arrays of the `feature_`/`algo_` columns to use. |
| `opts.sifted.flag` | `true` | Automatic feature selection (SIFTED). |
| `opts.sifted.K` | `10` | Number of features SIFTED selects. |
| `opts.pilot.dims` | `2` | Dimension of the instance space, 2 or 3. |
| `opts.pilot.method` | `'standard'` | Projection: `'standard'` or `'pls'` (Partial Least Squares). |
| `opts.pythia.classifier` | `'knn'` | Algorithm-selection classifier: `'knn'`, `'svm'`, `'tree'`, `'nb'`, `'linear'` or `'ensemble'`. |
| `opts.pythia.tuning` | `'sobol'` | Hyperparameter search: `'sobol'`, `'bayes'` or `'none'`. |
| `opts.trace.PI` | `0.6` | Minimum purity of a footprint. |
| `opts.outputs.csv`, `opts.outputs.png` | `true` | Write CSV files and PNG figures. |

The option groups follow the pipeline order: `general` and `perf` apply throughout; `prelim`, `auto`, `bound`, `norm` and `selvars` control data preparation (PRELIM); then `sifted` (feature selection), `pilot` (projection), `cloister` (boundary estimation), `pythia` (algorithm selection), `trace` (footprints), and `outputs` (files written at the end). An invalid value, such as `opts.pilot.dims = 4`, raises an `ISA:ISAvalidateOpts:*` error when the `InstanceSpace` object is created.

## AI-assisted analysis

This repository ships a [Claude Code](https://claude.com/claude-code) skill at `.claude/skills/instance-space-analysis/` to help an AI assistant perform Instance Space Analysis with this MATLAB toolkit — interpreting PRELIM/SIFTED/PILOT/CLOISTER/PYTHIA/TRACE output, choosing and justifying options (e.g. the good-performance threshold, classifier, feature-selection settings), designing a new ISA application, or debugging a pipeline run. It combines a general methodology reference (`SKILL.md`, covering the six-space framework and the six-step ISA process against the source papers) with an operational reference for this specific codebase (`references/matlab-toolkit.md`, covering the repository layout, `options.json` schema, and the `InstanceSpace` class API, kept current against this repository's actual behaviour rather than reconstructed from the papers alone).

**This skill only documents the MATLAB toolkit in this repository** — it does not cover MATILDA's web interface or the Python `pyInstanceSpace`/`matilda` package, whose options and code details differ from what's described here. The general methodology reference (six-space framework, reading footprints, when to augment vs. stop) is not toolkit-specific and may still help interpret results produced by those other tools, but any option name, file format, or code-behaviour detail should be verified against the tool actually in use.

## Contact

If you have any suggestions or ideas (e.g. for new features), or if you encounter any problems while running the code, please use the [issue tracker](https://github.com/andremun/InstanceSpace/issues) or contact us through MATILDA's [Queries and Feedback](http://matilda.unimelb.edu.au/matilda/contact-us) page. See [`CONTRIBUTING.md`](CONTRIBUTING.md) before opening a pull request, and [`SECURITY.md`](SECURITY.md) to report a suspected vulnerability privately.

## Acknowledgements

Funding for the development of this code was provided by:

- The Australian Research Council, through the Australian Laureate Fellowship FL140100012.
- The Australian Research Council, through the ARC Industrial Transformation Training Centre in Optimisation Technologies, Integrated Methodologies, and Applications (OPTIMA); grant No. IC200100009.
- The University of Melbourne, through grant 2025DYA013.
