# Instance Space Analysis: A toolkit for the assessment of algorithmic power

[![View InstanceSpace on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://au.mathworks.com/matlabcentral/fileexchange/75170-instancespace)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4750845.svg)](https://doi.org/10.5281/zenodo.4750845)

Instance Space Analysis is a methodology for assessing the strengths and weaknesses of an algorithm and objectively comparing its algorithmic power, without bias introduced by a restricted choice of test instances. At its core is the modelling of the relationship between an instance's structural properties and the performance of a group of algorithms. Instance Space Analysis allows the construction of **footprints** for each algorithm, defined as regions in the instance space where we statistically infer good performance. Other insights that can be gathered from Instance Space Analysis include:

-	Objective metrics of each algorithm’s footprint across the instance space as a measure of algorithmic power;
-	Explanation through visualisation of how instance features correlate with algorithm performance in various regions of the instance space;
-	Visualisation of the distribution and diversity of existing benchmark and real-world instances;
-	Assessment of the adequacy of the features used to characterise an instance;
-	Partitioning of the instance space into recommended regions for automated algorithm selection;
-	Distinguishing areas of the instance space where it may be useful to generate additional instances to gain further insights.

The unique advantage of visualising algorithm performance in the instance space, rather than as a small set of summary statistics averaged across a selected collection of instances, is the nuanced analysis that becomes possible to explain strengths and weaknesses and examine interesting variations in performance that may be hidden by tables of summary statistics.

This repository provides a set of MATLAB tools for performing a complete Instance Space Analysis as part of an automated pipeline. It is also the computational engine that powers the Melbourne Algorithm Test Instance Library with Data Analytics ([MATILDA](http://matilda.unimelb.edu.au/matilda/)) web tools for online analysis. For further information on the Instance Space Analysis methodology, see [here](http://matilda.unimelb.edu.au/matilda/our-methodology).

If you follow the Instance Space Analysis methodology, please cite as follows:

> K. Smith-Miles and M.A. Muñoz. *Instance Space Analysis for Algorithm Testing: Methodology and Software Tools*. ACM Comput. Surv. 55(12:255),1-31 [DOI:10.1145/3572895](https://doi.org/10.1145/3572895), 2023.

If you use the 3D extension of the methodology (```opts.pilot.dims = 3```, ```PILOTviewpoint```, TRACE3's native 3D footprints), please additionally cite:

> C. Simpson, M.A. Muñoz, S. Kandanaarachchi and R.J.G.B. Campello. *ISA3: A 3-dimensional expansion of Instance Space Analysis*. Machine Learning, 114, 240 [DOI:10.1007/s10994-025-06871-5](https://doi.org/10.1007/s10994-025-06871-5), 2025.

Also, if you specifically use this code, please cite as follows:

> M.A. Muñoz and K. Smith-Miles. *Instance Space Analysis: A toolkit for the assessment of algorithmic power*. andremun/InstanceSpace on GitHub. Zenodo, [DOI:10.5281/zenodo.4484107](https://doi.org/10.5281/zenodo.4484107), 2020.

Or if you specifically use [MATILDA](http://matilda.unimelb.edu.au/matilda/), please cite as follows:

> K. Smith-Miles, M.A. Muñoz and Neelofar. *Melbourne Algorithm Test Instance Library with Data Analytics (MATILDA)*. Available at (https://matilda.unimelb.edu.au). 2020.

**DISCLAIMER: This repository contains research code. On occasion, new features will be added, or changes will be made that may result in crashes. Although we have made every effort to reduce bugs, this code comes with NO GUARANTEES. If you find any issues, let us know ASAP via the contact methods listed at the end of this document.**

## Installation Instructions

The main requirement for the software to run is to have MATLAB R2025a or later, with the [Global Optimization](https://au.mathworks.com/help/gads/index.html), [Parallel Computing](https://www.mathworks.com/products/parallel-computing.html), [Optimization](https://au.mathworks.com/products/optimization.html), and [Statistics and Machine Learning](https://au.mathworks.com/help/stats/index.html) toolboxes installed. The Communications and Financial toolboxes are **not** required (an earlier version used the Communications Toolbox's `de2bi`; it has since been replaced with a pure-MATLAB equivalent). LIBSVM support is deprecated: new runs always use MATLAB's native classifier registry (```opts.pythia.classifier```). The LIBSVM MEX-files (`svmpredict`/`svmtrain`) are **not bundled with this repository** — they were precompiled binaries with no corresponding source in the tree, so they were removed (see #29). They're only relevant at all for evaluating a pre-v1.7 model whose classifiers `ISAmigrateModel` couldn't retrain (missing original training data); `ISAmigrateModel` prefers retraining from scratch with the current registry whenever the training data is available, which needs no LIBSVM dependency. If you do hit that path, obtain LIBSVM yourself from [the official project](https://www.csie.ntu.edu.tw/~cjlin/libsvm/) and add its MEX-files to the MATLAB path — `PYTHIA`'s eval mode raises a clear `ISA:PYTHIA:noLibsvm` error naming the missing dependency if you don't.

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
```

```InstanceSpace.m```/```buildIS.m```/```exploreIS.m``` add ```core/```, ```output/```, ```utils/```, and ```deprecated/``` to the MATLAB path automatically the first time any of them is used in a session — ```example.m```, ```test_integration.m```, and any script that starts with ```buildIS```/```exploreIS```/```InstanceSpace``` need no extra setup. If you want to call a function from one of those folders directly (```PILOT```, ```SIFTED```, ...) without going through one of those three first, run ```startup.m``` (e.g. ```run('startup.m')```, or just ```startup``` with the repo root as your current folder) at the start of your session.

## Working with the code

Start with ```example.m```: it runs the full pipeline (```buildIS``` + ```exploreIS```) once, on the bundled reference dataset, with sensible defaults and just a handful of commonly adjusted settings (classifier, tuning strategy, projection dimensionality, feature selection on/off) exposed as plain variables near the top. Outputs — images (```.png```), tables (```.csv```), and raw intermediate data (```.mat```) — land in ```test/data/example/```. To analyse your own data, point it at a folder containing your ```metadata.csv``` instead (see "The metadata file" below), and revisit the performance-metric settings, which are tuned for the bundled dataset's error-rate semantics.

```test_integration.m``` is the exhaustive option-coverage regression suite used during development — every classifier, tuning strategy, 2D/3D, PLS, viewpoint groups, and more, each in its own subdirectory under ```test/data/``` (e.g. ```test/data/classifier_svm/```) so no run overwrites another's outputs. It's a good reference for how a given option is meant to be used, but not the place to start.

**```options.json``` is a generated artifact, not a source file**, for both scripts above. Each run writes its own ```options.json``` from the ```opts``` struct built in MATLAB. Hand-editing an ```options.json``` file has no lasting effect — the next run silently overwrites it. To change what gets run, edit the MATLAB script instead (```example.m``` directly, or for ```test_integration.m```, the ```defaultOpts()``` local function for shared settings and a specific test case's ```override``` function for that case only).

### The InstanceSpace class

```buildIS```/```exploreIS``` are thin backward-compatibility wrappers (kept for callers like the MATILDA web platform that invoke them directly) around the ```InstanceSpace``` class, which is the recommended interface for new code:

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

See the class's own help text (```help InstanceSpace```) for the full method list, including ```plot()``` and ```getResults(idx)``` for accessing a specific ```explore()``` call's results.

## The metadata file

The ```metadata.csv``` file should contain a table where each row corresponds to a problem instance, and each column must strictly follow the naming convention mentioned below:

-	**instances** instance identifier - We expect the instance identifier to be of type "String". This column is mandatory.
-	**source** instance source - This column is optional
-	**feature_name** The keyword "feature_" concatenated with feature name. For instance, if the feature name is "density", the header name should be mentioned as "feature_density". If the name consists of more than one word, each word should be separated by "_" (spaces are not allowed). There must be more than two features for the software to work. We expect the features to be of the type "Double".
-	**algo_name** The keyword "algo_" concatenated with algorithm name. For instance, if the algorithm name is "Greedy", the column header should be "algo_greedy". If the name consists of more than one word, each word should be separated by "_" (spaces are not allowed). You can add the performance of more than one algorithm in the same ```.csv```. We expect the algorithm performance to be of the type "Double".

Moreover, empty cells, NaN or null values are allowed but **not recommended**. We expect you to handle missing values in your data before processing. You may use [this file](https://matilda.unimelb.edu.au/matilda/matildadata/graph_coloring_problem/metadata/metadata.csv) as reference.

**Common data-preparation mistake**: using `NA` instead of `NaN`, or leaving Excel error codes (`#REF!`, `#NULL!`, `#DIV/0!`) or empty rows in the sheet. Any of these causes `readtable` to infer a column as text (`string`/`cell`) instead of numeric (`double`), which will crash the pipeline downstream rather than failing with a clear error at load time.

## Options

Every setting below is a field of the ```opts``` structure passed to ```buildIS``` (as ```options.json```, built by ```example.m``` or ```test_integration.m```). The sections follow the pipeline's actual execution order (```InstanceSpace.build()```'s canonical stage order and dependencies, spec §7.3): general settings apply throughout; PRELIM (data bounding/scaling) runs first, then SIFTED (feature selection), then PILOT (dimensionality reduction); CLOISTER (bound estimation) and PYTHIA (algorithm selection) both run next, since each depends only on PILOT's output, not on each other; TRACE (footprint construction) runs last among the analysis stages, since it depends on PYTHIA's predictions; output settings apply to what's written once the pipeline completes.

### General settings

-	```opts.general.seed``` master RNG seed (default ```42```). Governs every stochastic stage of the pipeline (PILOT's BFGS random starts, PYTHIA's Sobol/Bayes tuning, SIFTED's GA, etc.), so a run with the same seed and inputs is reproducible.
-	```opts.general.verbose``` turns on (```TRUE```, default) or off (```FALSE```) the detailed, stage-level console output (per-trial/per-iteration progress, hyperparameter results, projection matrix display). One-line stage start/complete messages are always printed regardless of this setting. ```opts.pilot.verbose``` and ```opts.pythia.verbose``` inherit this value by default and can be overridden independently if needed.
- ```opts.general.parallel``` determines whether parallel processing will be available (set as ```TRUE```), or not (set as ```FALSE```). The toolkit uses MATLAB's [```parpool```](https://au.mathworks.com/help/parallel-computing/parpool.html) functionality to create a multi-session environment on the local machine.
- ```opts.general.ncores``` number of available cores for parallel processing.
-	```opts.perf.MaxPerf``` determines whether the algorithm performance values provided are **efficiency** measures that should be maximised (set as ```TRUE```), or **cost** measures that should be minimised (set as ```FALSE```).
-	```opts.perf.AbsPerf``` determines whether good performance is defined absolutely, e.g., misclassification error is lower than a 20%, (set as ```TRUE```), or if it is defined relatively to the best performing algorithm, e.g., misclassification error is within at least 5% of the best algorithm, (set as ```FALSE```).
-	```opts.perf.epsilon``` corresponds to the threshold used to calculate good performance. It must be of the type "Double".
-	```opts.perf.betaThreshold``` corresponds to the fraction of algorithms in the portfolio that must have good performance in the instance, for it to be considered an **easy** instance. It must be a value between 0 and 1.
-	```opts.selvars.feats``` / ```opts.selvars.algos``` optional cell arrays of feature/algorithm names to restrict the analysis to a hand-picked subset, matching the ```metadata.csv``` column headers **with** their ```feature_```/```algo_``` prefix (e.g. ```{'feature_density', 'feature_diameter'}```, not ```{'density', 'diameter'}```). Omit either field to use all available features/algorithms.
-	```opts.selvars.smallscaleflag``` By setting this flag as ```TRUE```, you can carry out a small-scale experiment using a randomly selected fraction of the original data. This is useful if you have a large dataset with more than 1000 instances and want to explore the model's parameters.
-	```opts.selvars.smallscale``` fraction taken from the original data on the small scale experiment.
-	```opts.selvars.fileidxflag``` by setting this flag as ```TRUE```, you can carry out a small scale experiment. This time you must provide a ```.csv``` file that contains in one column the indices of the instances to be taken. This may be useful if you want to make a more controlled experiment than just randomly selecting instances.
-	```opts.selvars.fileidx``` name of the file containing the indexes of the instances.
-	```opts.selvars.densityflag``` by setting this flag as ```TRUE```, instances are subset by feature-space density via ```FILTER``` instead of the small-scale/file-index options above: pairs of instances closer than ```opts.selvars.mindistance``` in feature space are treated as redundant, and one of each such pair is dropped, keeping a representative spread rather than a uniform random sample.
-	```opts.selvars.mindistance``` feature-space distance threshold below which two instances are considered too close (redundant) for the density-based filter.
-	```opts.selvars.type``` selects which extra condition, on top of feature-space closeness, ```FILTER``` requires before treating an instance as redundant and dropping it: ```'Ftr'``` (feature-space closeness alone), ```'Ftr&AP'``` (also requires similar algorithm performance), ```'Ftr&Good'``` (default; also requires both instances to be uniformly good across the whole portfolio), or ```'Ftr&AP&Good'``` (all of the above).

### Automatic data bounding and scaling

The toolkit implements simple routines to bound outliers and scale the data. **These routines are by no means perfect, and users should pre-process their data independently if preferred**. However, the automatic bounding and scaling routines should give some idea of the kind of results may be achieved. In general, we recommend transforming the data to be **close to normally distributed** due to the linear nature of PILOT's optimal projection algorithm.

- ```opts.auto.preproc``` turns on (set as ```TRUE```) the automatic pre-processing.
- ```opts.bound.flag``` turns on (set as ```TRUE```) data bounding. This sub-routine calculates the median and the interquartile range ([IQR](https://en.wikipedia.org/wiki/Interquartile_range)) of each feature and performance measure, and bounds the data to the median plus or minus ```opts.prelim.iqrMultiplier``` (default 5) times the IQR. **Warning**: this sub-routine often has issues with features that have low value diversity (e.g. most instances share the exact same value), since the IQR can collapse to zero. Consider either removing such features or setting ```opts.bound.flag = false```.
- ```opts.norm.flag``` turns on (set as ```TRUE```) scaling. This sub-routine scales each feature and performance measure into a positive range. Then it calculates a [box-cox transformation](https://en.wikipedia.org/wiki/Power_transform#Box%E2%80%93Cox_transformation) to stabilise the variance, and a [Z-transformation](https://en.wikipedia.org/wiki/Standard_score) to standardise the data. The results are features and performance measures that are close to normally distributed.
- ```opts.prelim.nanThreshold``` (default 0.20) maximum fraction of missing (```NaN```) values allowed for a feature before it is dropped entirely.

### Automatic feature selection

The toolkit implements SIFTED (```SIFTED.m```; ```SIFTED2``` is a deprecated alias kept for backward compatibility), a routine for selecting features based on their cross-correlation and correlation with performance. Ideally, we want the fewest orthogonal and predictive features. **This routine is by no means perfect, and users should pre-process their data independently if preferred**. In general, we recommend **using no more than 10 features** as input to PILOT's optimal projection algorithm, due to the numerical nature of its solution and issues in identifying meaningful linear trends.

- ```opts.sifted.flag``` turns on (set as ```TRUE```) the automatic feature selection. SIFTED runs in two stages:
	1. **Correlation filter.** SIFTED calculates the [Pearson correlation coefficient](https://en.wikipedia.org/wiki/Pearson_correlation_coefficient) between each feature and each algorithm's performance, keeping the single most-correlated feature for every algorithm unconditionally, plus any other feature whose absolute correlation with any algorithm exceeds ```opts.sifted.rho``` at significance ```opts.sifted.pval```. If 3 or fewer features remain at this point (or arrive with 3 or fewer to begin with), SIFTED stops here.
	2. **Correlation clustering + GA.** The Pearson correlation coefficient is used as a dissimilarity metric between the surviving features, and [k-means clustering](https://en.wikipedia.org/wiki/K-means_clustering) (```opts.sifted.K``` clusters, ```opts.sifted.MaxIter```/```opts.sifted.Replicates``` controlling convergence) groups similar features together. A [Genetic Algorithm](https://en.wikipedia.org/wiki/Genetic_algorithm) then searches for the best one-feature-per-cluster combination: each candidate combination is scored by projecting it with PILOT's analytic branch (a single closed-form solve, no restarts needed) at the same dimensionality as the outer pipeline's ```opts.pilot.dims```, then taking the worst-case (maximum) k-fold cross-validated KNN classification loss across all algorithms — lower is better. Previously-evaluated combinations are cached for the duration of the call. This stage is skipped if 3 or fewer features survive the correlation filter, or if there are no more features than ```opts.sifted.K```.
- ```opts.sifted.rho``` correlation threshold indicating the lowest acceptable absolute correlation between a feature and performance. It should be a value between 0 and 1.
- ```opts.sifted.pval``` significance level for the correlation filter; a feature/algorithm correlation is only counted if its p-value is at or below this threshold.
- ```opts.sifted.K``` number of clusters, which corresponds to the final number of features returned. The routine assumes at least 3 clusters and no more than the number of features. Ideally, it **should not** exceed 10.
- ```opts.sifted.MaxIter``` number of iterations used to converge the k-means algorithm. Usually, this setting does not need tuning.
- ```opts.sifted.Replicates``` number of repeats carried out of the k-means algorithm. Usually, this setting does not need tuning.

### Dimensionality reduction settings

The toolkit uses PILOT for dimensionality reduction and [BFGS](https://en.wikipedia.org/wiki/Broyden-Fletcher-Goldfarb-Shanno_algorithm) as a numerical solver. Technical details about it can be found [here](https://doi.org/10.1007/s10994-017-5629-5).

-	```opts.pilot.analytic``` determines whether the analytic (set as ```TRUE```) or the numerical (set as ```FALSE```) solution to the dimensionality reduction problem should be used. We recommend leaving this setting as ```FALSE```due to the instability of the analytical solution caused by possible poor conditioning. Only applies when ```opts.pilot.method = 'standard'```.
-	```opts.pilot.ntries``` number of iterations that the numerical solution is attempted.
-	```opts.pilot.dims``` projection dimensionality, ```2``` (default) or ```3```. Replaces the legacy boolean ```opts.pilot.ISA3D```, which is still accepted as an alias (```ISA3D = true``` maps to ```dims = 3```). A 3D projection lets you additionally call ```PILOTviewpoint``` to find the best 2D camera angle(s) onto it (see below).
-	```opts.pilot.method``` selects the projection algorithm: ```'standard'``` (default) is the BFGS/analytic method described above; ```'pls'``` uses Partial Least Squares ([```plsregress```](https://au.mathworks.com/help/stats/plsregress.html)) instead, which maximises covariance between the projection and the performance matrix and does not require the feature matrix to be full column rank. Both methods work at 2D or 3D via ```opts.pilot.dims```.
-	```opts.pilot.alpha``` (default ```1.0```) scales the performance-reconstruction term of PILOT's cost function relative to the feature-reconstruction term, for ```opts.pilot.method = 'standard'``` only: `min ||F̃-BrZ||² + α||Y-CrZ||²`. Increase it to emphasise performance trends over feature trends in the projection.
-	```opts.pilot.topoWeight``` (default ```0```, disabled) is reserved for future experimental use and has no effect in the current version.

When ```opts.pilot.dims = 3```, ```buildIS``` automatically calls ```PILOTviewpoint(Z, Y, opts.pilot)``` to find the best 2D camera viewpoint(s) of the 3D projection, storing the result in ```model.pilot.viewpoint```. By default, one viewpoint is found across all algorithms; ```opts.pilot.viewGroups``` (a cell array of algorithm index vectors, e.g. ```{[1 2 3], [4 5 6]}```) requests one viewpoint per group instead, useful for inspecting a subset of algorithms in isolation.

### Empirical bound estimation settings

The toolkit uses CLOISTER, an algorithm based on correlation to detect the empirical bounds of the Instance Space. It runs after PILOT, on the projected instance space, and does not depend on PYTHIA.

- ```opts.cloister.corrThreshold``` Determines the maximum [Pearson correlation coefficient](https://en.wikipedia.org/wiki/Pearson_correlation_coefficient) that would indicate non-correlated variables. The lower this value is, the more stringent the algorithm is; hence, it would be less likely to produce a good bound.
- ```opts.cloister.pval``` Determines the p-value of the Pearson correlation coefficient that indicates no correlation.
- ```opts.cloister.maxFeatures``` (default ```20```) guards the number of features CLOISTER's correlation-based bit-matrix enumeration will process; above this count it would become intractable, so CLOISTER instead falls back to a plain convex hull of the projected instances as the boundary, with a warning.

### Algorithm selection settings

The toolkit selects one binary classifier per algorithm (good/not-good performance) from a registry of MATLAB-native classifiers, resolved via `ISAgetClassifierFcn`. It runs after PILOT, on the projected instance space, and does not depend on CLOISTER. **LIBSVM is deprecated for new runs**: `buildIS` never dispatches to it. `ISAmigrateModel` (see below) renames legacy field names on an old model (e.g. `.svm`/`.knn` → `.classifiers`), and additionally **retrains** any classifiers still in the legacy LIBSVM struct format (no `predict()` method, so they can't just be relabelled) using the current registry — `opts.pythia.classifier` if already set to a valid registry name, `'knn'` otherwise — provided the model still has `model.pilot.Z` and `model.data.Y`/`Ybin`/`Ybest`/`algolabels` to retrain from. If those fields are missing, migration proceeds but leaves the legacy classifiers in place with a warning; such a model can still be evaluated through `exploreIS`/`PYTHIA` eval mode (which dispatches to `svmpredict` when it detects a struct instead of a MATLAB classifier object), but this requires the LIBSVM MEX-files, which are **not bundled with this repository** (see the Installation Instructions section above) — obtain them from [the official LIBSVM project](https://www.csie.ntu.edu.tw/~cjlin/libsvm/) if you actually need to go down this path. Prefer retraining from scratch with `buildIS`, or a migration with the full training data available, to drop that dependency entirely.

- ```opts.pythia.classifier``` selects the classifier: ```'knn'``` (default, via `fitcknn`), ```'svm'``` (`fitcsvm`), ```'tree'``` (`fitctree`), ```'nb'``` (Naive Bayes, `fitcnb`), ```'linear'``` (`fitclinear`), or ```'ensemble'``` (`fitcensemble`; see ```opts.pythia.ensembleMethod```, default ```'Bag'```). All algorithms in the portfolio use the same classifier.
- ```opts.pythia.tuning``` selects the hyperparameter search strategy: ```'sobol'``` (default; a scrambled Sobol quasi-random sequence, ```opts.pythia.nTuningIter``` evaluations with k-fold CV), ```'bayes'``` (MATLAB `bayesopt`, Gaussian-process surrogate, same evaluation budget and CV), or ```'none'``` (use ```opts.pythia.params``` directly, skipping tuning).
- ```opts.pythia.nTuningIter``` number of Sobol/Bayes evaluations (default 20).
- ```opts.pythia.kFold``` number of folds of the CV experiment.
- ```opts.pythia.params``` pre-calculated hyperparameters; required when ```opts.pythia.tuning = 'none'```, and always takes precedence over tuning when supplied. Shape depends on ```opts.pythia.classifier```: ```[nalgos x 1]``` for the single-parameter classifiers (```'tree'```, ```'nb'```, ```'linear'```), or ```[nalgos x 2]``` for the two-parameter ones (```'knn'```, ```'svm'```, ```'ensemble'```).
- ```opts.pythia.skip``` bypasses classifier training entirely (TRACE then falls back to the true labels directly, with a warning).
- ```opts.pythia.ispolykrnl``` (SVM only) determines whether to use a polynomial (set as ```TRUE```) or Gaussian (set as ```FALSE```) kernel. Usually, the latter is significantly faster to compute and more accurate; however, it also has the disadvantage of producing discontinuous areas of good performance, which may appear overfit. We tend to recommend a polynomial kernel if the dataset is larger than 1000 instances.
- ```opts.pythia.useweights``` determines whether weighted (set as ```TRUE```) or unweighted (set as ```FALSE```) classification is performed. The weights are calculated as <img src="https://render.githubusercontent.com/render/math?math=\left|y_{\text{best}}-y\right|">.
- ```opts.pythia.seed``` RNG seed for classifier training/tuning; defaults to ```opts.general.seed``` and rarely needs to be set independently.
- ```opts.pythia.verbose``` per-classifier training progress and tuning output; defaults to ```opts.general.verbose```.

### Footprint construction settings

The toolkit uses TRACE3, an algorithm based on MATLAB's [```alphaShape```](https://au.mathworks.com/help/matlab/ref/alphashape.html) to define the regions in the space where we statistically infer good algorithm performance, applicable to both 2D and 3D instance spaces. It runs last, after PYTHIA: TRACE3 always reuses PYTHIA's predicted labels for the good-performance region (`Zu = {z_i : yhat_i=1 AND ybin_i=1}`) rather than retraining its own classifier — this coupling is unconditional and not configurable. When ```opts.pythia.skip = true```, TRACE falls back to the true labels only (`Zu = {z_i : ybin_i=1}`), with a warning.

- ```opts.trace.method``` selects ```'trace3'``` (default, above) or ```'legacy'``` (the pre-refactor DBSCAN + alpha-shape triangulation method, 2D only).
- ```opts.trace.PI``` minimum purity required for a section of a footprint.
- ```opts.trace.minInstances``` (default ```4```) minimum number of instances a candidate footprint must contain to be considered valid.
- ```opts.trace.minAreaFrac``` (default ```0.01```) minimum footprint size, as a fraction of the total instance-space area/volume, for a candidate footprint to be kept.
- ```opts.trace.contra``` (```'legacy'``` method only, defaults to ```TRUE``` there) turns on contradiction removal — trimming overlapping sections of two algorithms' footprints where the evidence is weak. Not read by the default ```'trace3'``` method.

### Output settings

These settings result in more information being stored in files or presented in the console output once the pipeline above has completed.

- ```opts.outputs.csv``` This flag produces the output CSV files for post-processing and analysis. It is recommended to leave this setting as ```TRUE```.
- ```opts.outputs.png``` This flag produces the output figure files for post-processing and analysis. It is recommended to leave this setting as ```TRUE```.
- ```opts.outputs.fig``` For 3D projections, this flag also writes a ```.fig``` file alongside each footprint PNG for interactive rotation in MATLAB (see ```ISArecallView``` to snap a reopened figure back to its optimised viewpoint). No effect for 2D projections.
- ```opts.outputs.web``` This flag produces the output files employed to draw the figures in MATILDA's web tools (click [here](https://matilda.unimelb.edu.au/matilda/newuser) to open an account). It is recommended to leave this setting as ```FALSE```.

## AI-assisted analysis

This repository ships a [Claude Code](https://claude.com/claude-code) skill at ```.claude/skills/instance-space-analysis/``` to help an AI assistant perform Instance Space Analysis with this MATLAB toolkit — interpreting PRELIM/SIFTED/PILOT/CLOISTER/PYTHIA/TRACE output, choosing and justifying options (e.g. the good-performance threshold, classifier, feature-selection settings), designing a new ISA application, or debugging a pipeline run. It combines a general methodology reference (```SKILL.md```, covering the six-space framework and the six-step ISA process against the source papers) with an operational reference for this specific codebase (```references/matlab-toolkit.md```, covering the repository layout, ```options.json``` schema, and the ```InstanceSpace``` class API, kept current against this repository's actual behaviour rather than reconstructed from the papers alone).

**This skill only documents the MATLAB toolkit in this repository** — it does not cover MATILDA's web interface or the Python ```pyInstanceSpace```/```matilda``` package, whose options and code details differ from what's described here. The general methodology reference (six-space framework, reading footprints, when to augment vs. stop) is not toolkit-specific and may still help interpret results produced by those other tools, but any option name, file format, or code-behaviour detail should be verified against the tool actually in use.

## Contact

If you have any suggestions or ideas (e.g. for new features), or if you encounter any problems while running the code, please use the [issue tracker](https://github.com/andremun/InstanceSpace/issues) or contact us through MATILDA's [Queries and Feedback](http://matilda.unimelb.edu.au/matilda/contact-us) page. See [`CONTRIBUTING.md`](CONTRIBUTING.md) before opening a pull request, and [`SECURITY.md`](SECURITY.md) to report a suspected vulnerability privately.

## Acknowledgements

Funding for the development of this code was provided by:

- The Australian Research Council, through the Australian Laureate Fellowship FL140100012.
- The Australian Research Council, through the ARC Industrial Transformation Training Centre in Optimisation Technologies, Integrated Methodologies, and Applications (OPTIMA); grant No. IC200100009.
- The University of Melbourne, through grant 2025DYA013.
