% example.m - Getting-started example for the ISA Toolkit.
%
% Runs the full Instance Space Analysis pipeline (buildIS + exploreIS) on
% the bundled reference dataset (Munoz et al. 2018 classification study:
% 212 instances x 10 features x 10 classification algorithms, scored by
% misclassification error) using sensible defaults, with just a handful
% of the settings most people adjust first exposed as plain variables
% below.
%
% To analyse your own data instead: replace SRCDIR with the folder
% containing your metadata.csv (see README.md, "The metadata file"), and
% revisit the "performance metric" settings below -- they are tuned for
% this dataset's error-rate semantics, not a generic default.
%
% For the exhaustive option-coverage regression suite (every classifier,
% tuning strategy, 2D/3D, PLS, viewpoint groups, etc.), see
% test_integration.m instead. For the full list of configurable options,
% see README.md, "Options".

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
%
% Reference:
%   Smith-Miles, K. & Munoz, M.A. (2023). Instance Space Analysis for
%   Algorithm Testing. ACM Computing Surveys, 55(12), Article 255.
%   https://doi.org/10.1145/3572895
% -------------------------------------------------------------------------

srcdir  = './test/data/';
rootdir = './test/data/example/';
if ~isfolder(rootdir), mkdir(rootdir); end
copyfile([srcdir 'metadata.csv'], [rootdir 'metadata.csv']);
copyfile([srcdir 'metadata_test.csv'], [rootdir 'metadata_test.csv']);

% ---- A few of the most commonly adjusted settings --------------------------
classifier = 'knn';    % opts.pythia.classifier: 'knn' (default), 'svm', 'tree', 'nb', 'linear', 'ensemble'
tuning     = 'sobol';  % opts.pythia.tuning: 'sobol' (default), 'bayes', or 'none' (needs opts.pythia.params)
dims       = 2;        % opts.pilot.dims: 2 (default) or 3
siftedFlag = true;     % opts.sifted.flag: automated feature selection on/off
% -----------------------------------------------------------------------

% Starts from a clean struct so a stale opts variable left in the base
% workspace by a previous script (e.g. test_integration.m) in the same
% MATLAB session can't leak unrelated fields into this run's options.json.
opts = struct();
opts.pythia.classifier = classifier;
opts.pythia.tuning     = tuning;
opts.pilot.dims        = dims;
opts.sifted.flag       = siftedFlag;

% This dataset's algorithm performance is a misclassification error rate:
% lower is better (MaxPerf=false), and "good" means an absolute error
% below 20% (AbsPerf=true, epsilon=0.20). Everything else is left at the
% toolkit's defaults (see ISAdefaults.m / README.md, "Options").
opts.perf.MaxPerf = false;
opts.perf.AbsPerf = true;
opts.perf.epsilon = 0.20;

fid = fopen([rootdir 'options.json'], 'w+');
fprintf(fid, '%s', jsonencode(opts));
fclose(fid);

fprintf('[EXAMPLE] Building the instance space...\n');
model = buildIS(rootdir); %#ok<NASGU>

fprintf('[EXAMPLE] Evaluating the held-out test instances...\n');
out = exploreIS(rootdir); %#ok<NASGU>

fprintf('[EXAMPLE] Done. Outputs (CSVs, PNGs, model.mat) are in %s\n', rootdir);
fprintf('[EXAMPLE] ''model'' and ''out'' are in your workspace to explore.\n');
fprintf('[EXAMPLE] See README.md for the full list of configurable options.\n');
