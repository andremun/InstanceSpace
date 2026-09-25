classdef ValidationTest < matlab.unittest.TestCase
% ValidationTest  Unit tests for the option and utility functions in
% utils/: every ISAvalidateOpts check, ISAdefaults' legacy-name mapping,
% the ISAgetClassifierFcn registry, and ISAsubsetData. These run on
% synthetic inputs in milliseconds, so they cover error branches that a
% full pipeline run never reaches.

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

    properties (TestParameter)
        BadOption = badOptionCases();
        Classifier = {'knn', 'svm', 'tree', 'nb', 'linear', 'ensemble'};
    end

    methods (TestClassSetup)
        function addToolkitPath(~)
            % These tests call utils/ functions directly, without first
            % creating an InstanceSpace object (which would add the path).
            root = testRepoRoot();
            addpath([root 'core'], [root 'utils'], [root 'output']);
        end
    end

    methods (Test)
        function testRejectsInvalidOption(testCase, BadOption)
            testCase.verifyError(@() ISAvalidateOpts(BadOption.opts), BadOption.id);
        end

        function testRejectsNonStruct(testCase)
            testCase.verifyError(@() ISAvalidateOpts(5), 'ISA:ISAvalidateOpts:notStruct');
        end

        function testAcceptsCompleteDefaults(testCase)
            % Every default must itself pass validation, or InstanceSpace.load
            % (which validates a stored, defaults-filled opts) would fail.
            opts = ISAdefaults(struct());
            testCase.verifyEqual(ISAvalidateOpts(opts), opts);
        end

        function testAcceptsValidEdgeValues(testCase)
            opts = struct();
            opts.general = struct('seed', 0, 'ncores', 'auto');     % 0 seed; non-numeric ncores
            opts.perf = struct('AbsPerf', true, 'epsilon', 25);      % absolute epsilon outside [0,1]
            opts.pilot = struct('dims', 3, 'topoWeight', 0, 'viewGroups', {{[1 2], 3}});
            opts.pythia = struct('classifier', "SVM", 'ensembleMethod', "Bag");
            opts.selvars = struct('feats', {{'feature_a', "feature_b"}}, 'fileidx', "idx.csv");
            testCase.verifyEqual(ISAvalidateOpts(opts), opts);
        end

        function testDefaultsFillEveryGroup(testCase)
            opts = ISAdefaults(struct());
            groups = {'general', 'perf', 'prelim', 'auto', 'bound', 'norm', 'selvars', ...
                      'sifted', 'pilot', 'cloister', 'pythia', 'trace', 'outputs'};
            for i = 1:numel(groups)
                testCase.verifyTrue(isstruct(opts.(groups{i})), ...
                    sprintf('ISAdefaults should create opts.%s.', groups{i}));
            end
            testCase.verifyEqual(opts.sifted.seed, opts.general.seed);
            testCase.verifyEqual(opts.pythia.seed, opts.general.seed);
            testCase.verifyFalse(opts.trace.contra, ...
                'Contradiction removal defaults to false for trace3.');
        end

        function testDefaultsKeepUserValues(testCase)
            opts.pilot.dims = 3;
            opts.general.seed = 7;
            opts = ISAdefaults(opts);
            testCase.verifyEqual(opts.pilot.dims, 3);
            testCase.verifyEqual(opts.pilot.seed, 7, 'Stage seeds inherit a user-set general seed.');
        end

        function testDefaultsMapLegacyNames(testCase)
            opts.parallel = struct('flag', true, 'ncores', 4);
            opts.pilot.ISA3D = true;
            opts.cloister.cthres = 0.8;
            opts.pythia.cvfolds = 7;
            opts.pythia.useknn = false;
            opts.trace.method = 'legacy';
            opts = ISAdefaults(opts);
            testCase.verifyTrue(opts.general.parallel);
            testCase.verifyEqual(opts.general.ncores, 4);
            testCase.verifyEqual(opts.pilot.dims, 3);
            testCase.verifyEqual(opts.cloister.corrThreshold, 0.8);
            testCase.verifyEqual(opts.pythia.kFold, 7);
            testCase.verifyEqual(opts.pythia.classifier, 'svm');
            testCase.verifyTrue(opts.trace.contra, ...
                'Contradiction removal defaults to true for the legacy method.');
        end

        function testClassifierRegistry(testCase, Classifier)
            [fitFcn, p1, p2] = ISAgetClassifierFcn(Classifier);
            testCase.verifyClass(fitFcn, 'function_handle');
            testCase.verifyTrue(startsWith(func2str(fitFcn), 'fitc'));
            testCase.verifyNotEmpty(p1);
            testCase.verifyNotEmpty(p2);
            % Names are case insensitive.
            testCase.verifyEqual(func2str(ISAgetClassifierFcn(upper(Classifier))), func2str(fitFcn));
        end

        function testClassifierRegistryRejectsUnknown(testCase)
            testCase.verifyError(@() ISAgetClassifierFcn('forest'), ...
                'ISA:ISAgetClassifierFcn:unknownClassifier');
        end

        function testMigrateModelFile(testCase)
            % File form of ISAmigrateModel: backup, write-back, and its guards.
            tmp = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            testCase.verifyError(@() ISAmigrateModel(tmp), 'ISA:ISAmigrateModel:noModelFile');

            legacy = struct('opts', struct('pbldr', struct('dims', 2), ...
                                           'perf', struct('MaxMin', true)));
            save(fullfile(tmp, 'model.mat'), '-struct', 'legacy');
            migrated = ISAmigrateModel(tmp);
            testCase.verifyTrue(isfield(migrated.opts, 'pilot') && ~isfield(migrated.opts, 'pbldr'));
            testCase.verifyTrue(migrated.opts.perf.MaxPerf);
            testCase.verifyTrue(isfile(fullfile(tmp, 'model_legacy.mat')), 'The original is backed up.');
            onDisk = load(fullfile(tmp, 'model.mat'));
            testCase.verifyTrue(isfield(onDisk.opts, 'pilot'), 'The migrated model is written back.');

            testCase.verifyError(@() ISAmigrateModel(tmp), 'ISA:ISAmigrateModel:backupExists', ...
                'An existing backup must never be overwritten.');
            ISAmigrateModel(tmp, 'backupSuffix', "_second");
            testCase.verifyTrue(isfile(fullfile(tmp, 'model_second.mat')));

            testCase.verifyWarning(@() ISAmigrateModel(legacy, 'backupSuffix', '_x'), ...
                'ISA:ISAmigrateModel:ignoredArgs');
        end

        function testSubsetData(testCase)
            data = syntheticData(6, 4, 3);
            data.S = categorical({'a'; 'b'; 'a'; 'b'; 'a'; 'b'});
            keep = logical([1 0 1 0 1 0])';
            out = ISAsubsetData(data, keep);
            testCase.verifySize(out.X, [3 4]);
            testCase.verifyEqual(out.instlabels, data.instlabels(keep));
            testCase.verifyEqual(out.S, data.S(keep));
            testCase.verifyEqual(out.P, data.P(keep));

            out = ISAsubsetData(data, [2 4], [1 3]);
            testCase.verifyEqual(out.X, data.X([2 4], [1 3]));
            testCase.verifyEqual(out.featlabels, data.featlabels([1 3]));
            testCase.verifyEqual(out.Xraw, data.Xraw([2 4], :), ...
                'Only X and featlabels are column-subset; Xraw keeps every feature.');
        end
    end
end

% =========================================================================
function cases = badOptionCases()
% One invalid value per check in ISAvalidateOpts, with the error it raises.
c = {
    'generalNotStruct',    'general',  [],              5,          'notStruct'
    'verboseNumeric',      'general',  'verbose',       1,          'notLogical'
    'parallelEmpty',       'general',  'parallel',      [],         'notLogical'
    'seedNegative',        'general',  'seed',          -1,         'notPositive'
    'seedFraction',        'general',  'seed',          1.5,        'notInteger'
    'ncoresZero',          'general',  'ncores',        0,          'notPositive'
    'maxPerfText',         'perf',     'MaxPerf',       'yes',      'notLogical'
    'epsilonRelative',     'perf',     'epsilon',       2,          'notInUnitRange'
    'betaNaN',             'perf',     'betaThreshold', NaN,        'notInUnitRange'
    'iqrZero',             'prelim',   'iqrMultiplier', 0,          'notPositive'
    'nanThresholdHigh',    'prelim',   'nanThreshold',  1.2,        'notInUnitRange'
    'preprocNumeric',      'auto',     'preproc',       0,          'notLogical'
    'boundText',           'bound',    'flag',          'on',       'notLogical'
    'normNumeric',         'norm',     'flag',          1,          'notLogical'
    'smallscaleHigh',      'selvars',  'smallscale',    2,          'notInUnitRange'
    'fileidxNumeric',      'selvars',  'fileidx',       5,          'notText'
    'mindistanceNegative', 'selvars',  'mindistance',   -0.1,       'notPositive'
    'filterType',          'selvars',  'type',          'Near',     'notMember'
    'featsNumeric',        'selvars',  'feats',         {1, 2},     'notCellOfText'
    'algosText',           'selvars',  'algos',         'algo_a',   'notCellOfText'
    'siftedRho',           'sifted',   'rho',           -0.5,       'notInUnitRange'
    'siftedK',             'sifted',   'K',             2.5,        'notInteger'
    'siftedReplicates',    'sifted',   'Replicates',    0,          'notPositive'
    'pilotMethod',         'pilot',    'method',        'pca',      'notMember'
    'pilotDims',           'pilot',    'dims',          4,          'notMember'
    'pilotDimsText',       'pilot',    'dims',          '3',        'notMember'
    'pilotNtries',         'pilot',    'ntries',        0,          'notPositive'
    'pilotAlphaZero',      'pilot',    'alpha',         0,          'notPositive'
    'pilotTopoNegative',   'pilot',    'topoWeight',    -1,         'notPositive'
    'viewGroupsNumeric',   'pilot',    'viewGroups',    [1 2],      'badViewGroups'
    'viewGroupsZero',      'pilot',    'viewGroups',    {[1 0]},    'badViewGroups'
    'cloisterPval',        'cloister', 'pval',          1.5,        'notInUnitRange'
    'cloisterMaxFeatures', 'cloister', 'maxFeatures',   -3,         'notPositive'
    'pythiaClassifier',    'pythia',   'classifier',    'forest',   'notMember'
    'pythiaTuning',        'pythia',   'tuning',        'grid',     'notMember'
    'pythiaIter',          'pythia',   'nTuningIter',   Inf,        'notPositive'
    'pythiaKFold',         'pythia',   'kFold',         2.5,        'notInteger'
    'pythiaEnsemble',      'pythia',   'ensembleMethod', 3,         'notText'
    'traceMethod',         'trace',    'method',        'dbscan',   'notMember'
    'tracePI',             'trace',    'PI',            1.1,        'notInUnitRange'
    'traceMinInstances',   'trace',    'minInstances',  0,          'notPositive'
    'traceContra',         'trace',    'contra',        'no',       'notLogical'
    'outputsCsv',          'outputs',  'csv',           [],         'notLogical'
    'outputsFig',          'outputs',  'fig',           2,          'notLogical'
    };
cases = struct();
for i = 1:size(c, 1)
    opts = struct();
    if isempty(c{i,3})
        opts.(c{i,2}) = c{i,4};
    else
        opts.(c{i,2}).(c{i,3}) = c{i,4};
    end
    cases.(c{i,1}) = struct('opts', opts, 'id', ['ISA:ISAvalidateOpts:' c{i,5}]);
end
% epsilon under AbsPerf=true must still be a finite real scalar.
cases.epsilonAbsoluteNaN = struct( ...
    'opts', struct('perf', struct('AbsPerf', true, 'epsilon', NaN)), ...
    'id', 'ISA:ISAvalidateOpts:notFiniteNumericScalar');
end

% =========================================================================
function data = syntheticData(ninst, nfeats, nalgos)
% Minimal model.data-shaped struct for ISAsubsetData.
rng(1, 'twister');
data.X = rand(ninst, nfeats);
data.Y = rand(ninst, nalgos);
data.Xraw = data.X + 1;
data.Yraw = data.Y + 1;
data.Ybin = data.Y > 0.5;
data.beta = any(data.Ybin, 2);
data.numGoodAlgos = sum(data.Ybin, 2);
data.Ybest = min(data.Y, [], 2);
[~, data.P] = min(data.Y, [], 2);
data.instlabels = arrayfun(@(i) sprintf('inst%d', i), (1:ninst)', 'UniformOutput', false);
data.featlabels = arrayfun(@(i) sprintf('f%d', i), 1:nfeats, 'UniformOutput', false);
end
