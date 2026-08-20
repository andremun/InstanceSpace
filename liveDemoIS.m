%% Instance Space Analysis: A toolkit for the assessment of algorithmic power
% liveDemoIS.m - Interactive, stage-by-stage walkthrough of the ISA
% pipeline using the |InstanceSpace| class.
%
% Instance Space Analysis is a methodology for the assessment of the
% strengths and weaknesses of an algorithm, and an approach to
% objectively compare algorithmic power without the bias introduced by a
% restricted choice of test instances. At its core is the modelling of
% the relationship between structural properties of an instance and the
% performance of a group of algorithms. Instance Space Analysis allows
% the construction of _footprints_ for each algorithm: regions in the
% instance space where we statistically infer good performance. Other
% insights include:
%
% * Objective metrics of each algorithm's footprint as a measure of
% algorithmic power;
% * Explanation, through visualisation, of how instance features
% correlate with algorithm performance in different regions of the
% instance space;
% * Visualisation of the distribution and diversity of existing
% benchmark and real-world instances;
% * Assessment of the adequacy of the features used to characterise an
% instance;
% * Partitioning of the instance space into recommended regions for
% automated algorithm selection;
% * Identification of areas of the instance space where generating
% additional instances would give the most insight.
%
% This script walks through the |InstanceSpace| class one stage at a
% time -- PRELIM, SIFTED, PILOT, CLOISTER, PYTHIA, then TRACE -- running
% each stage, inspecting what it produced, and (where useful) plotting
% it, before moving to the next. It replaces the previous |liveDemoIS.mlx|
% Live Script, which predates the |InstanceSpace| class and several
% since-completed pipeline changes (PRELIM as its own function, PYTHIA's
% classifier registry, TRACE3, PILOT's 3D/PLS support, and more). It's
% shipped as a plain, |%%|-sectioned |.m| file -- MATLAB's Live Editor
% opens it exactly like a Live Script (Home > Open > this file, or drag
% it into an open Live Editor tab) and you can "Save As... > MATLAB Live
% Code File" from there if you specifically want an |.mlx| -- because
% plain text diffs and reviews far better in git than the Live Script's
% binary format does.
%
% *Not* every |opts| field is demonstrated here; only the handful with
% the most visible effect on the result. For the complete option
% reference, see README.md, "Options". For the exhaustive
% option-coverage regression suite (every classifier, tuning strategy,
% 2D/3D, PLS, viewpoint groups, etc.), see test_integration.m instead --
% this script is meant to teach the pipeline, not exercise every corner
% of it.
%
% If you use Instance Space Analysis, please cite:
%
% K. Smith-Miles & M.A. Munoz (2023). Instance Space Analysis for
% Algorithm Testing: Methodology and Software Tools. ACM Computing
% Surveys, 55(12), Article 255. https://doi.org/10.1145/3572895
%
% If you use this toolkit itself, please also cite:
%
% M.A. Munoz & K. Smith-Miles. Instance Space Analysis: A toolkit for
% the assessment of algorithmic power. andremun/InstanceSpace on GitHub.
% Zenodo, https://doi.org/10.5281/zenodo.4484107.
%
% *DISCLAIMER*: this toolkit is in active development. We've made every
% effort to keep it correct and well-tested (see test_integration.m),
% but it comes with no guarantees. If you find an issue, please open one
% on GitHub.

% -------------------------------------------------------------------------
% Instance Space Analysis (ISA) Toolkit
% Copyright (c) 2026 Mario Andres Munoz Acosta and contributors
% School of Computing and Information Systems
% The University of Melbourne, Australia
%
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
% License: https://polyformproject.org/licenses/noncommercial/1.0.0/
%
% You may use, modify, and redistribute this software for non-commercial
% research and educational purposes only. Commercial use requires prior
% written permission. See the LICENSE file for full terms.
% -------------------------------------------------------------------------

%% Installation instructions
% You need MATLAB R2025a or later, with the
% <https://au.mathworks.com/help/gads/index.html Global Optimization>,
% <https://www.mathworks.com/products/parallel-computing.html Parallel Computing>,
% <https://au.mathworks.com/products/optimization.html Optimization>,
% <https://au.mathworks.com/help/stats/index.html Statistics and Machine
% Learning>, and <https://au.mathworks.com/products/finance.html Financial>
% toolboxes installed (Financial Toolbox is needed for |boxcox()|, used by
% |PRELIM.m|'s auto-normalisation step). The Communications Toolbox is
% *not* required. LIBSVM support is deprecated: new runs
% always use MATLAB's native classifier registry
% (|opts.pythia.classifier|). The LIBSVM MEX-files are not bundled with
% this repository; they're only relevant for evaluating a model migrated
% from a pre-refactor toolkit version whose classifiers couldn't be
% retrained (see |ISAmigrateModel|), and can be obtained from the official
% LIBSVM project if actually needed.
%
% Run |startup.m| once per session before calling a |core|/|output|/|utils|
% function directly; the class constructor (used below) adds those
% folders to the path automatically the first time it runs.

%% Setting up the environment and loading the training data
% The |InstanceSpace| class takes a |rootdir| containing |metadata.csv|
% (see README.md, "The metadata file") and an optional |opts| struct
% (any field you don't set falls back to a documented default -- see
% |ISAdefaults.m| / README.md, "Options"). We use the bundled reference
% dataset: the classification study of Munoz et al. (2018, _Machine
% Learning_, 107(1), 109-147) -- 212 training instances, 23 held-out test
% instances, 10 features, 10 classification algorithms scored by
% misclassification error.
srcdir  = './test/data/';
rootdir = './test/data/liveDemo/';
if ~isfolder(rootdir), mkdir(rootdir); end
copyfile([srcdir 'metadata.csv'], [rootdir 'metadata.csv']);
copyfile([srcdir 'metadata_test.csv'], [rootdir 'metadata_test.csv']);

opts = struct();
% This dataset's algorithm performance is a misclassification error
% rate: lower is better (MaxPerf=false), and "good" means an absolute
% error below 20% (AbsPerf=true, epsilon=0.20).
opts.perf.MaxPerf = false;
opts.perf.AbsPerf = true;
opts.perf.epsilon = 0.20;

obj = InstanceSpace(rootdir, opts);
disp(obj.completedStages);   % {} -- nothing has run yet; the constructor
                              % only reads and validates metadata.csv/opts.

%% PRELIM: loading and preparing the data
% PRELIM removes instances/features with too many missing values,
% optionally bounds outliers (median +/- |opts.prelim.iqrMultiplier| x
% IQR) and normalises (Box-Cox + Z-score) every feature and performance
% column, then works out which algorithm is "good" on which instance
% (|opts.perf.*| above) -- this is what the rest of the pipeline builds
% on.
obj = obj.build('stages', {'prelim'});
disp(obj.completedStages);
fprintf('Instances: %d, features: %d, algorithms: %d\n', ...
    size(obj.model.data.X,1), size(obj.model.data.X,2), size(obj.model.data.Y,2));
fprintf('Fraction of instances that are "easy" (a good algorithm exists): %.2f\n', ...
    mean(obj.model.data.beta));

%% SIFTED: automated feature selection
% SIFTED keeps the features most predictive of algorithm performance --
% a correlation filter first, then (if more than a handful of features
% survive) a genetic algorithm that picks one representative feature per
% correlation cluster. Set |opts.sifted.flag = false| to skip it and use
% every feature instead.
obj = obj.build('stages', {'sifted'});
fprintf('Features before SIFTED: %s\n', strjoin(obj.model.featsel.labels, ', '));
fprintf('Features after  SIFTED: %s\n', strjoin(obj.model.data.featlabels, ', '));

%% PILOT: obtaining a two-dimensional projection
% PILOT finds a linear projection |Z = X*A''| of the (post-SIFTED)
% feature matrix that jointly reconstructs the features and the
% algorithm performance as closely as possible -- this |Z| is the
% instance space itself. |opts.pilot.dims = 3| would give a 3D
% projection instead (see |PILOTviewpoint| for finding a good camera
% angle onto it); |opts.pilot.method = 'pls'| swaps in Partial Least
% Squares as an alternative to the default BFGS/analytic method.
obj = obj.build('stages', {'pilot'});
figure;
scatter(obj.model.pilot.Z(:,1), obj.model.pilot.Z(:,2), 20, ...
    obj.model.data.numGoodAlgos, 'filled');
xlabel('z_1'); ylabel('z_2'); colorbar;
title('Instance space, coloured by number of good algorithms per instance');

%% CLOISTER: finding the empirical bounds of the space
% CLOISTER estimates the boundary of the region of the instance space
% that's actually reachable given the observed feature correlations --
% useful for judging whether a new instance (or a synthetically
% generated one) falls inside or outside the range the model was built
% from.
obj = obj.build('stages', {'cloister'});
figure;
obj.plot('boundary');

%% PYTHIA: building an oracle for algorithm selection
% PYTHIA trains one binary classifier per algorithm (good/not-good
% performance) over the instance space |Z|, using the registry selected
% by |opts.pythia.classifier| (|'knn'| by default; also |'svm'|,
% |'tree'|, |'nb'|, |'linear'|, |'ensemble'|), tuned via
% |opts.pythia.tuning| (a scrambled Sobol sequence by default; also
% |'bayes'| or pre-supplied |opts.pythia.params|).
obj = obj.build('stages', {'pythia'});
disp(obj.model.pythia.summary);   % accuracy/precision/recall per algorithm

figure;
scatter(obj.model.pilot.Z(:,1), obj.model.pilot.Z(:,2), 20, ...
    obj.model.pythia.selection0, 'filled');
xlabel('z_1'); ylabel('z_2'); colorbar;
title('PYTHIA''s predicted best algorithm per instance');

%% TRACE: calculating the algorithm footprints
% TRACE3 builds a footprint -- a region of the instance space where an
% algorithm is expected to perform well -- from PYTHIA's predictions,
% for both individual algorithms and the portfolio as a whole. This is
% the step |opts.trace.method = 'legacy'| would swap for the
% pre-refactor DBSCAN + alpha-shape algorithm instead (2D only).
obj = obj.build('stages', {'trace'});
disp(obj.model.trace.summary);   % area/density/purity per algorithm

%% Post-processing: preparing results for further analysis or publication
% |obj.save()| writes |model.mat| (HDF5-compatible, |-v7.3|) to
% |rootdir|, containing every stage's full output -- including the
% fitted classifiers and the footprint geometry objects -- for
% downstream analysis, plotting, or containment queries
% (|isinterior|/|inShape|). CSV and PNG outputs (|opts.outputs.csv|/|.png|,
% both on by default) are written as a side effect of |build()| itself,
% not |save()|.
obj.save();
fprintf('Outputs written to %s (model.mat, CSVs, PNGs).\n', rootdir);

%% Exploring the Instance Space model
% |explore()| projects new instances into the already-fitted space and
% evaluates PYTHIA/TRACE on them, without retraining anything. It reads
% |metadata_test.csv| from the directory you pass it -- here, the same
% |rootdir| the training data came from, which is why we copied
% |metadata_test.csv| there at the start.

%% Setting up the test data, and testing PYTHIA and the footprints together
obj = obj.explore(rootdir);
testResults = obj.getResults(1);   % first explore() call's results
disp(testResults.pythia.summary);
disp(testResults.trace.summary);

figure;
scatter(testResults.pilot.Z(:,1), testResults.pilot.Z(:,2), 30, ...
    testResults.pythia.selection0, 'filled');
xlabel('z_1'); ylabel('z_2'); colorbar;
title('Held-out test instances projected into the trained instance space');

%% Additional resources
% * README.md -- the complete option reference and repository layout.
% * test_integration.m -- exhaustive option-coverage regression suite;
% a good reference for how any specific option is meant to be used.
% * <https://matilda.unimelb.edu.au MATILDA> -- the web platform this
% toolkit's |buildIS|/|exploreIS| entry points feed into.
% * K. Smith-Miles, M.A. Munoz & Neelofar. Melbourne Algorithm Test
% Instance Library with Data Analytics (MATILDA). Available at
% https://matilda.unimelb.edu.au.

%% Acknowledgements
% This toolkit and the methodology it implements are the product of
% many years of collaborative research by the MATILDA team,
% The University of Melbourne, and collaborators. See README.md,
% "Acknowledgements" for the full list.
