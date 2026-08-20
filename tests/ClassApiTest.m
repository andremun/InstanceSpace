classdef ClassApiTest < matlab.unittest.TestCase
% ClassApiTest  InstanceSpace class API coverage (spec Phase 7): staged
% build() with an option change between stages, out-of-order stage
% requests, stage invalidation on re-running an earlier stage, the
% missing-prerequisite and not-built error paths, a save()/load()
% round-trip, and explore() without going through exploreIS.m.
%
% Migrated from test_integration.m's 'class_api' block (#39). Kept as one
% sequential test method (not split into several independent test*
% methods): each step depends on obj's state left by the previous one, so
% splitting would either need matlab.unittest's execution-order guarantees
% (which the framework deliberately doesn't provide) or a shared fixture
% rebuilt at every step -- neither is simpler than one method here.

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
        CaseDir
        BaseOpts
    end

    methods (TestClassSetup)
        function setupCaseDir(testCase)
            % Absolute, not './test/data/': matlab.unittest does not
            % guarantee the working directory during test execution is
            % the directory test_integration.m was launched from (see
            % testRepoRoot.m).
            rootdir = [testRepoRoot() 'test/data/'];
            srcfiles = {'metadata.csv', 'metadata_test.csv'};
            testCase.CaseDir = [rootdir 'class_api/'];
            if ~isfolder(testCase.CaseDir), mkdir(testCase.CaseDir); end
            for f = 1:numel(srcfiles)
                copyfile([rootdir srcfiles{f}], [testCase.CaseDir srcfiles{f}]);
            end
            testCase.BaseOpts = testDefaultOpts();
        end
    end

    methods (Test)
        function testStagedBuildSaveLoadExplore(testCase)
            classCaseDir = testCase.CaseDir;
            baseOpts = testCase.BaseOpts;

            obj = InstanceSpace(classCaseDir, baseOpts);

            obj = obj.build('stages', {'prelim', 'sifted', 'pilot'});
            testCase.verifyEqual(obj.completedStages, {'prelim', 'sifted', 'pilot'}, ...
                'completedStages mismatch after a partial build().');

            obj.opts.pilot.alpha = 2.0;
            obj = obj.build('stages', {'pilot'}); % re-run just PILOT with the new weight
            testCase.verifyEqual(obj.completedStages, {'prelim', 'sifted', 'pilot'}, ...
                're-running an already-completed stage should not duplicate it in completedStages.');

            % Requested out of canonical order: build() must still run them
            % prelim->...->trace internally regardless of the order listed here.
            obj = obj.build('stages', {'trace', 'cloister', 'pythia'});
            testCase.verifyTrue(all(ismember({'cloister', 'pythia', 'trace'}, obj.completedStages)), ...
                'cloister/pythia/trace did not complete.');

            % Re-running an EARLIER stage after later ones have already completed
            % must invalidate cloister/pythia/trace, not leave them looking valid
            % alongside a freshly re-run pilot.
            obj.opts.pilot.alpha = 3.0;
            obj = obj.build('stages', {'pilot'});
            testCase.verifyEqual(obj.completedStages, {'prelim', 'sifted', 'pilot'}, ...
                're-running pilot after cloister/pythia/trace completed should invalidate them in completedStages.');
            testCase.verifyFalse(isfield(obj.model, 'cloist') || isfield(obj.model, 'pythia') || isfield(obj.model, 'trace'), ...
                're-running pilot should remove the now-stale cloister/pythia/trace model fields.');

            % explore() must refuse a model left partially invalidated like this,
            % rather than crash deep inside evaluateTestSet on a missing field.
            testCase.verifyError(@() obj.explore(classCaseDir), 'ISA:InstanceSpace:notBuilt', ...
                'explore() on a partially-invalidated model should raise ISA:InstanceSpace:notBuilt.');

            % Complete the pipeline again before the save()/load()/explore() checks below.
            obj = obj.build('stages', {'cloister', 'pythia', 'trace'});

            % Re-running 'cloister' must NOT invalidate 'pythia'/'trace': per
            % StagePrereq neither depends on cloister's output (both depend on
            % 'pilot'/'pythia' instead), even though both appear later than
            % 'cloister' in canonical StageOrder.
            pythiaBefore = obj.model.pythia;
            traceBefore = obj.model.trace;
            obj = obj.build('stages', {'cloister'});
            testCase.verifyTrue(all(ismember({'cloister', 'pythia', 'trace'}, obj.completedStages)), ...
                're-running cloister should not invalidate pythia/trace in completedStages.');
            testCase.verifyTrue(isequal(obj.model.pythia, pythiaBefore) && isequal(obj.model.trace, traceBefore), ...
                're-running cloister should not recompute/discard the still-valid pythia/trace results.');

            % A genuinely missing prerequisite must error clearly, not crash deep
            % inside the requested stage.
            testCase.verifyError(@() InstanceSpace(classCaseDir, baseOpts).build('stages', {'pilot'}), ...
                'ISA:InstanceSpace:missingPrereq', ...
                'build(''stages'',{''pilot''}) on a fresh object should raise ISA:InstanceSpace:missingPrereq.');

            % save()/load() round-trip.
            obj.save();
            loaded = InstanceSpace.load(classCaseDir);
            testCase.verifyEqual(loaded.model.pilot.A, obj.model.pilot.A, ...
                'save()/load() round-trip changed model.pilot.A.');

            % explore() directly through the class, not via exploreIS.m.
            loaded = loaded.explore(classCaseDir);
            testCase.verifyTrue(numel(loaded.testResults) == 1 && strcmp(loaded.testDirs{1}, classCaseDir), ...
                'explore() did not record testResults/testDirs correctly.');
            testCase.verifyTrue(isfield(loaded.getResults(1), 'trace'), ...
                'getResults(1) is missing the trace field.');
        end
    end
end

