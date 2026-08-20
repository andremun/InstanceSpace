classdef MigrationTest < matlab.unittest.TestCase
% MigrationTest  ISAmigrateModel legacy-migration-table coverage (spec
% §6.4), plus the LIBSVM-removal (#29) and Yraw (#42) regressions that
% ride on the same trained model this class already needs to build for
% the LIBSVM-retraining/legacy-TRACE-recompute cases.
%
% Migrated from test_integration.m's 'isamigratemodel' block (#39). The
% pure opts/data-struct renames (no trained model needed) are independent
% test* methods; TestClassSetup builds one real trained model, shared via
% BaseModel, for the cases that need it -- replacing the pre-#39 script's
% "reuse obj from the class_api case if available, else build one"
% workaround with a proper shared fixture.

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

    properties
        BaseModel
    end

    methods (TestClassSetup)
        function buildBaseModel(testCase)
            % Absolute, not './test/data/': matlab.unittest does not
            % guarantee the working directory during test execution is
            % the directory test_integration.m was launched from (see
            % testRepoRoot.m).
            rootdir = [testRepoRoot() 'test/data/'];
            srcfiles = {'metadata.csv', 'metadata_test.csv'};
            caseDir = [rootdir 'migration_test/'];
            if ~isfolder(caseDir), mkdir(caseDir); end
            for f = 1:numel(srcfiles)
                copyfile([rootdir srcfiles{f}], [caseDir srcfiles{f}]);
            end
            migrateObj = InstanceSpace(caseDir, testDefaultOpts()).build();
            testCase.BaseModel = migrateObj.model;
        end
    end

    methods (Test)
        function testOptsStructRenames(testCase)
            m = ISAmigrateModel(struct('opts', struct( ...
                'pbldr',     struct('dims', 2), ...
                'sbound',    struct('pval', 0.05), ...
                'footprint', struct('method', 'trace3'))));
            testCase.verifyTrue(isfield(m.opts, 'pilot') && m.opts.pilot.dims == 2 && ~isfield(m.opts, 'pbldr'), ...
                'opts.pbldr -> opts.pilot rename failed.');
            testCase.verifyTrue(isfield(m.opts, 'cloister') && ~isfield(m.opts, 'sbound'), ...
                'opts.sbound -> opts.cloister rename failed.');
            testCase.verifyTrue(isfield(m.opts, 'trace') && ~isfield(m.opts, 'footprint'), ...
                'opts.footprint -> opts.trace rename failed.');
        end

        function testCorrClustMerge(testCase)
            m = ISAmigrateModel(struct('opts', struct( ...
                'corr',  struct('flag', true, 'threshold', 0.35), ...
                'clust', struct('flag', false))));
            testCase.verifyTrue(~isfield(m.opts, 'corr') && ~isfield(m.opts, 'clust'), ...
                'opts.corr/opts.clust were not removed after merging into opts.sifted.');
            testCase.verifyTrue(isfield(m.opts, 'sifted') && m.opts.sifted.rho == 0.35, ...
                'opts.corr.threshold did not map to opts.sifted.rho.');
            testCase.verifyEqual(m.opts.sifted.flag, false, ...
                'opts.clust.flag=false did not set opts.sifted.flag=false.');
        end

        function testMaxMinRename(testCase)
            m = ISAmigrateModel(struct('opts', struct('perf', struct('MaxMin', true))));
            testCase.verifyTrue(isfield(m.opts.perf, 'MaxPerf') && m.opts.perf.MaxPerf == true && ...
                ~isfield(m.opts.perf, 'MaxMin'), 'opts.perf.MaxMin -> opts.perf.MaxPerf rename failed.');
        end

        function testBestPerformaceTypoRename(testCase)
            m = ISAmigrateModel(struct('data', struct('bestPerformace', [1;2;3])));
            testCase.verifyTrue(isfield(m.data, 'Ybest') && isequal(m.data.Ybest, [1;2;3]) && ...
                ~isfield(m.data, 'bestPerformace'), 'model.data.bestPerformace -> model.data.Ybest rename failed.');
        end

        function testCompletedStagesInferred(testCase)
            m = ISAmigrateModel(struct('prelim', struct(), 'sifted', struct(), 'pilot', struct()));
            testCase.verifyEqual(m.completedStages, {'prelim','sifted','pilot'}, ...
                'completedStages was not correctly inferred from present model sub-structs.');
        end

        function testIncompletePilotWarning(testCase)
            testCase.verifyWarning(@() ISAmigrateModel(struct('pilot', struct('A', 1))), ...
                'ISA:ISAmigrateModel:incompletePilot', ...
                'model.pilot.A without B/C should raise an ISA:ISAmigrateModel:incompletePilot warning.');
        end

        function testLibsvmRetrainingAndYrawRegression(testCase)
            baseModel = testCase.BaseModel;
            % LIBSVM struct: a plain struct (no predict() method) under the old
            % 'svm' field name, as svmtrain() would have produced. NOTE: assigned
            % directly (pythiaLegacy.svm = {...}), not via struct('svm',{...}) --
            % the latter's cell-value form builds a struct ARRAY from the cell's
            % elements instead of a scalar struct whose field holds the cell array.
            pythiaLegacy = struct();
            pythiaLegacy.svm = repmat({struct('nr_class', 2)}, 1, numel(baseModel.data.algolabels));
            libsvmModel = baseModel;
            libsvmModel.pythia = pythiaLegacy;
            m = ISAmigrateModel(libsvmModel);
            testCase.verifyTrue(isfield(m.pythia, 'classifiers') && ~isfield(m.pythia, 'svm'), ...
                'LIBSVM-format pythia.svm was not retrained/migrated into pythia.classifiers.');
            testCase.verifyEqual(numel(m.pythia.classifiers), numel(baseModel.data.algolabels), ...
                'Retrained pythia.classifiers has the wrong number of algorithms.');

            % #42 regression: retrainLibsvmPythia previously passed model.data.Y
            % (normalized) instead of model.data.Yraw to PYTHIA.
            % Avg/Std_Perf_all_instances are pure functions of the Y passed in
            % (independent of which classifier was trained), so a
            % normalized-scale bug here would show up as a scale mismatch
            % against baseModel's own summary, built through the normal path
            % on the same data.
            avgPerfRetrained = cell2mat(m.pythia.summary(2:end-2, 2));
            avgPerfOriginal  = cell2mat(baseModel.pythia.summary(2:end-2, 2));
            testCase.verifyEqual(avgPerfRetrained, avgPerfOriginal, ...
                ['Retrained pythia.summary''s Avg_Perf_all_instances does not match a freshly-built ' ...
                 'model''s summary on the same data -- retrainLibsvmPythia is training on the wrong ' ...
                 '(normalized) scale.']);
        end

        function testNoLibsvmRegression(testCase)
            % #29 regression: the LIBSVM MEX-files (svmpredict/svmtrain) were
            % removed from this repository as precompiled binaries with no
            % build source in the tree. PYTHIA's eval mode must now give a
            % clear, actionable ISA:PYTHIA:noLibsvm error instead of MATLAB's
            % generic "undefined function 'svmpredict'" if it ever needs to
            % dispatch to a legacy (un-retrained) LIBSVM-format classifier
            % struct. Assumes no LIBSVM MEX-files are on this MATLAB path,
            % which holds for this repository post-#29 unless one is
            % installed separately.
            baseModel = testCase.BaseModel;
            pythiaLegacy = struct();
            pythiaLegacy.svm = repmat({struct('nr_class', 2)}, 1, numel(baseModel.data.algolabels));
            testCase.verifyError(@() PYTHIA(baseModel.pilot.Z, baseModel.data.Yraw, baseModel.data.Ybin, ...
                baseModel.data.Ybest, baseModel.data.algolabels, baseModel.opts.pythia, pythiaLegacy), ...
                'ISA:PYTHIA:noLibsvm', ...
                ['PYTHIA eval mode on a legacy LIBSVM struct should raise ISA:PYTHIA:noLibsvm when ' ...
                 'svmpredict is unavailable, not MATLAB''s generic undefined-function error.']);
        end

        function testLegacyTraceRecompute(testCase)
            % Legacy DBSCAN+polyshape triangulation format: .space.area, no .measure.
            baseModel = testCase.BaseModel;
            legacyTraceModel = baseModel;
            legacyTraceModel.trace = struct('space', struct('area', 1, 'polygon', []));
            m = ISAmigrateModel(legacyTraceModel);
            testCase.verifyTrue(isfield(m.trace.space, 'measure') && ~isfield(m.trace.space, 'area'), ...
                'Legacy triangulation-format model.trace was not recomputed via TRACE3.');
            testCase.verifyTrue(isfield(m.trace, 'good') && numel(m.trace.good) == numel(baseModel.data.algolabels), ...
                'Recomputed model.trace is missing per-algorithm footprints.');
        end
    end
end
