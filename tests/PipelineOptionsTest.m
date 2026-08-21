classdef PipelineOptionsTest < matlab.unittest.TestCase
% PipelineOptionsTest  Option-surface coverage: every classifier, tuning
% strategy, 2D/3D, PLS, viewpoint groups, etc., run through buildIS()+
% exploreIS() against the bundled reference dataset.
%
% Migrated from test_integration.m's testCases cell array (#39): each
% entry there is now a case of the OptionCase TestParameter below, run by
% the single parameterized test method testOption. A failing case reports
% its own name via the standard TestParameter diagnostic instead of a
% hand-rolled pass/fail struct array.
%
% IMPORTANT: options.json is NOT the source of truth for these settings.
% It is a generated artifact: TestClassSetup writes a fresh options.json
% into each case's subdirectory from the opts struct built in MATLAB
% below. Hand-editing an options.json file has no lasting effect -- the
% next run overwrites it. To change what gets tested, edit
% testDefaultOpts() (tests/testDefaultOpts.m) or the relevant OptionCase
% entry.

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

    properties (Constant)
        % Absolute, not './test/data/': matlab.unittest does not guarantee
        % the working directory during test execution is the directory
        % test_integration.m was launched from (see testRepoRoot.m).
        RootDir = [testRepoRoot() 'test/data/'];
        SrcFiles = {'metadata.csv', 'metadata_test.csv'};
    end

    properties (TestParameter)
        % Each field name is the case name (also the output subdirectory,
        % matching the pre-#39 convention so test/data/<case_name>/ output
        % paths are unchanged). Value is a struct with .desc (what the
        % case exercises) and .override (function handle applied to
        % testDefaultOpts()).
        OptionCase = pipelineOptionCases();
    end

    methods (TestClassSetup)
        function ensureRootDir(testCase) %#ok<MANU>
            if ~isfolder(PipelineOptionsTest.RootDir)
                mkdir(PipelineOptionsTest.RootDir);
            end
        end
    end

    methods (Test)
        function testOption(testCase, OptionCase)
            fprintf('[TEST] %s\n', OptionCase.desc);
            caseDir = [PipelineOptionsTest.RootDir OptionCase.name '/'];
            if ~isfolder(caseDir), mkdir(caseDir); end
            for f = 1:numel(PipelineOptionsTest.SrcFiles)
                copyfile([PipelineOptionsTest.RootDir PipelineOptionsTest.SrcFiles{f}], ...
                    [caseDir PipelineOptionsTest.SrcFiles{f}]);
            end

            opts = OptionCase.override(testDefaultOpts());
            fid = fopen([caseDir 'options.json'], 'w+');
            fprintf(fid, '%s', jsonencode(opts));
            fclose(fid);

            % No try/catch, no verifyWarningFree: an uncaught error here
            % fails this test automatically, which is all the original
            % script's try/catch bookkeeping ever checked for. Several
            % cases (e.g. pythia_skip) legitimately raise a warning as
            % part of correct behaviour, so warnings must not be treated
            % as failures.
            runCase(caseDir);
        end
    end
end

function runCase(caseDir)
model = buildIS(caseDir); %#ok<NASGU>
out = exploreIS(caseDir); %#ok<NASGU>
end

function cases = pipelineOptionCases()
% Number of algorithms in the reference dataset, read from the CSV header
% rather than hardcoded, so the 'tuning_none_precalc' case's pre-supplied
% opts.pythia.params matrix is always sized correctly even if the dataset
% changes. TestParameter properties are evaluated once, statically, before
% any test runs, so the CSV must already exist in the repo (it does --
% test/data/metadata.csv is tracked, see the repo's .gitignore notes).
hdr = readtable([testRepoRoot() 'test/data/metadata.csv']);
nalgos = sum(strncmpi(hdr.Properties.VariableNames, 'algo_', 5));

cases = struct();

cases.baseline_knn_sobol_2d = caseEntry('baseline_knn_sobol_2d', ...
    ['Default pipeline: KNN classifier, Sobol tuning, TRACE3, SIFTED on, 2D PILOT. ' ...
     'The happy path every other case is a variation of.'], ...
    @(o) o);

cases.classifier_svm = caseEntry('classifier_svm', ...
    'PYTHIA classifier registry: SVM (fitcsvm) instead of the KNN default.', ...
    @(o) setPythia(o, 'classifier', 'svm'));

cases.classifier_ensemble = caseEntry('classifier_ensemble', ...
    ['PYTHIA classifier registry: ensemble (fitcensemble, Bag method) -- the ' ...
     'two-hyperparameter non-KNN case (NumLearningCycles/MinLeafSize).'], ...
    @(o) setPythia(setPythia(o, 'classifier', 'ensemble'), 'ensembleMethod', 'Bag'));

cases.tuning_bayes = caseEntry('tuning_bayes', ...
    ['opts.pythia.tuning=''bayes'': MATLAB bayesopt hyperparameter search instead ' ...
     'of the Sobol default.'], ...
    @(o) setPythia(o, 'tuning', 'bayes'));

cases.tuning_none_precalc = caseEntry('tuning_none_precalc', ...
    ['opts.pythia.tuning=''none'' with pre-supplied opts.pythia.params: skips ' ...
     'hyperparameter search entirely and uses fixed values.'], ...
    @(o) setPythiaPrecalc(o, nalgos));

cases.pythia_skip = caseEntry('pythia_skip', ...
    ['opts.pythia.skip=true: PYTHIA bypasses classifier training; TRACE must fall ' ...
     'back to the true-label-only footprint (Zu={yi=1}) with a warning.'], ...
    @(o) setPythia(o, 'skip', true));

cases.trace_legacy = caseEntry('trace_legacy', ...
    'opts.trace.method=''legacy'': the pre-refactor DBSCAN + alpha-shape footprint algorithm.', ...
    @(o) setTrace(o, 'method', 'legacy'));

cases.sifted_off = caseEntry('sifted_off', ...
    'opts.sifted.flag=false: pipeline runs without automated feature selection.', ...
    @(o) setSifted(o, 'flag', false));

cases.small_scale_subset = caseEntry('small_scale_subset', ...
    'opts.selvars.smallscaleflag=true: validates the small-scale subsetting path.', ...
    @(o) setSmallScale(o, 0.5));

cases.density_subset = caseEntry('density_subset', ...
    'opts.selvars.densityflag=true: validates the FILTER-based density subsetting path.', ...
    @(o) setDensitySubset(o));

cases.cost_sensitive_weights = caseEntry('cost_sensitive_weights', ...
    ['opts.pythia.useweights=true: cost-sensitive classification weighted by ' ...
     '|ybest-y|, including the degenerate-weights fallback guard.'], ...
    @(o) setPythia(o, 'useweights', true));

cases.pilot_3d_numeric = caseEntry('pilot_3d_numeric', ...
    ['opts.pilot.dims=3, opts.pilot.analytic=false: the numerical (BFGS) 3D ' ...
     'projection path. Also exercises PILOTviewpoint''s default global-viewpoint ' ...
     'path (buildIS calls it automatically whenever dims=3) and scriptpng''s use ' ...
     'of the resulting azimuth/elevation for every 3D plot (opts.outputs.png=true ' ...
     'by default).'], ...
    @(o) setPilot(o, 'dims', 3));

cases.pilot_analytic_2d = caseEntry('pilot_analytic_2d', ...
    'opts.pilot.analytic=true, 2D: the closed-form eigenvector solution (top-2 eigenvectors).', ...
    @(o) setPilot(o, 'analytic', true));

cases.pilot_analytic_3d = caseEntry('pilot_analytic_3d', ...
    ['opts.pilot.analytic=true, opts.pilot.dims=3: the generalised closed-form ' ...
     'solution (top-3 eigenvectors) -- exercises the 3D analytic branch specifically.'], ...
    @(o) setPilot(setPilot(o, 'analytic', true), 'dims', 3));

cases.pilot_isa3d_legacy_alias = caseEntry('pilot_isa3d_legacy_alias', ...
    ['opts.pilot.ISA3D=true with dims removed: regression check that the legacy ' ...
     'boolean alias still maps to the 3D path via ISAdefaults/PILOT''s backward-compat shim.'], ...
    @(o) setPilotISA3DLegacy(o));

cases.pilot_pls_2d = caseEntry('pilot_pls_2d', ...
    'opts.pilot.method=''pls'', 2D: the Partial Least Squares alternative projection method.', ...
    @(o) setPilot(o, 'method', 'pls'));

cases.pilot_pls_3d = caseEntry('pilot_pls_3d', ...
    'opts.pilot.method=''pls'', opts.pilot.dims=3: PLS at 3D.', ...
    @(o) setPilot(setPilot(o, 'method', 'pls'), 'dims', 3));

cases.pilot_alpha_weight = caseEntry('pilot_alpha_weight', ...
    ['opts.pilot.alpha=3.0: emphasises the performance-reconstruction term over the ' ...
     'feature-reconstruction term in PILOT''s cost function.'], ...
    @(o) setPilot(o, 'alpha', 3.0));

cases.pilot_viewpoint_groups = caseEntry('pilot_viewpoint_groups', ...
    ['opts.pilot.dims=3, opts.pilot.viewGroups={1:4,5:10}: exercises PILOTviewpoint''s ' ...
     'per-group path (one viewpoint per algorithm group instead of a single global ' ...
     'one), and scriptpng''s per-algorithm group lookup (resolveViewAngle) that picks ' ...
     'the right viewpoint for each algorithm''s plots. Uneven group sizes so the ' ...
     'options.json round-trip decodes as a cell array rather than collapsing into a ' ...
     'numeric matrix (jsonencode/jsondecode quirk).'], ...
    @(o) setPilot(setPilot(o, 'dims', 3), 'viewGroups', {1:4, 5:10}));

cases.feats_algos_subset = caseEntry('feats_algos_subset', ...
    ['opts.selvars.feats/algos: restrict to a hand-picked subset of features and ' ...
     'algorithms (README-documented, must be preserved).'], ...
    @(o) setFeatsAlgosSubset(o));

cases.sifted_parallel_ga = caseEntry('sifted_parallel_ga', ...
    ['opts.general.parallel=true: exercises SIFTED''s GA-level UseParallel path and ' ...
     'clearCache''s cross-worker persistent-cache reset via parfevalOnAll (spec Phase 6). ' ...
     'Every other case leaves opts.general.parallel=false, so this path would ' ...
     'otherwise go untested.'], ...
    @(o) setGeneralParallel(o));
end

function s = caseEntry(name, desc, override)
% name is stored explicitly (not inferred from the containing struct's
% field name) so testOption can use it directly -- e.g. as the case's
% output subdirectory -- without depending on matlab.unittest's
% TestParameter-name reflection inside the test method itself.
s = struct('name', name, 'desc', desc, 'override', override);
end

function opts = setPythia(opts, field, value)
opts.pythia.(field) = value;
end

function opts = setSifted(opts, field, value)
opts.sifted.(field) = value;
end

function opts = setPilot(opts, field, value)
opts.pilot.(field) = value;
end

function opts = setGeneralParallel(opts)
opts.general.parallel = true;
opts.general.ncores = 2; % small pool, just enough to exercise the parallel path
end

function opts = setPilotISA3DLegacy(opts)
% Removes the dims field (set unconditionally in testDefaultOpts) so
% ISAdefaults' legacy-alias migration -- which only fires when dims is
% absent -- actually gets exercised by this case.
opts.pilot = rmfield(opts.pilot, 'dims');
opts.pilot.ISA3D = true;
end

function opts = setTrace(opts, field, value)
opts.trace.(field) = value;
end

function opts = setPythiaPrecalc(opts, nalgos)
opts.pythia.tuning = 'none';
% Fixed KNN hyperparameters for every algorithm: NumNeighbors=5, Distance
% index 1 (euclidean, per ISAgetClassifierFcn's distance ordering).
% nalgos is the raw algorithm count from metadata.csv's header; this
% assumes buildIS.m's "algorithms with no good instances" filter doesn't
% drop any of them under this case's perf/epsilon settings. If it ever
% does, PYTHIA's own validation raises a clear ISA:PYTHIA:noParamsForNoneTuning
% error rather than failing silently.
opts.pythia.params = [repmat(5, nalgos, 1), repmat(1, nalgos, 1)];
end

function opts = setSmallScale(opts, fraction)
opts.selvars.smallscaleflag = true;
opts.selvars.smallscale = fraction;
end

function opts = setDensitySubset(opts)
opts.selvars.densityflag = true;
opts.selvars.mindistance = 0.1;
opts.selvars.type = 'Ftr&Good';
end

function opts = setFeatsAlgosSubset(opts)
% Names must match the raw CSV header (feature_/algo_ prefixes intact):
% buildIS.m applies this filter before stripping the prefixes.
opts.selvars.feats = {'feature_Max_Normalized_Entropy_attributes', ...
                       'feature_ErrorRate_Decision_Node', ...
                       'feature_Max_Feature_Efficiency_F3', ...
                       'feature_Training_Error_Linear_Classifier_L2', ...
                       'feature_Nonlinearity_Nearest_Neighbor_Classifier_N4'};
opts.selvars.algos = {'algo_NB', 'algo_CART', 'algo_KNN'};
end
