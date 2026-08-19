classdef RegressionTest < matlab.unittest.TestCase
% RegressionTest  Targeted regressions for specific fixed bugs (#41, #44,
% #28, #37/#38), each verifying a fix rather than exercising an option.
% Migrated from test_integration.m's bespoke post-testCases blocks (#39).

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
        BaseModel
    end

    methods (TestClassSetup)
        function setupFixture(testCase)
            rootdir = './test/data/';
            srcfiles = {'metadata.csv', 'metadata_test.csv'};
            testCase.CaseDir = [rootdir 'regression_test/'];
            if ~isfolder(testCase.CaseDir), mkdir(testCase.CaseDir); end
            for f = 1:numel(srcfiles)
                copyfile([rootdir srcfiles{f}], [testCase.CaseDir srcfiles{f}]);
            end
            testCase.BaseOpts = testDefaultOpts();
            testCase.BaseModel = InstanceSpace(testCase.CaseDir, testCase.BaseOpts).build().model;
        end
    end

    methods (Test)
        function testSeedReproducibility(testCase)
            % #41: PILOT/SIFTED previously always reseeded via
            % rng('default'), ignoring opts.general.seed entirely -- a user
            % rerunning with a different seed for a replication/variance
            % study got identical PILOT/SIFTED results anyway. Confirms a
            % different seed changes PILOT's BFGS-multi-start result, and
            % that the same seed still reproduces bit-identical results.
            classCaseDir = testCase.CaseDir;
            baseOpts = testCase.BaseOpts;

            seedOptsA = baseOpts;
            seedOptsA.general.seed = 1;
            objSeedA = InstanceSpace(classCaseDir, seedOptsA).build('stages', {'prelim', 'sifted', 'pilot'});

            seedOptsB = baseOpts;
            seedOptsB.general.seed = 2;
            objSeedB = InstanceSpace(classCaseDir, seedOptsB).build('stages', {'prelim', 'sifted', 'pilot'});

            % PILOT's numerical (BFGS) path draws its multi-start X0 from
            % rng(opts.seed,'twister'); different seeds must give different
            % random starting points and therefore (with overwhelming
            % probability, given the non-convex cost surface) a different
            % converged A. This is the robust half of the check -- SIFTED's
            % K-means+GA feature selection can legitimately converge to the
            % same selvars across seeds on a small dataset even when its own
            % seeding is correct, so that's logged for information rather
            % than asserted.
            testCase.verifyNotEqual(objSeedA.model.pilot.A, objSeedB.model.pilot.A, ...
                'Different opts.general.seed values produced identical PILOT.A -- the seed is not threading through PILOT''s BFGS multi-start.');
            if isequal(objSeedA.model.sifted.selvars, objSeedB.model.sifted.selvars)
                fprintf('[TEST] Note: SIFTED.selvars happened to match across seeds 1 and 2 (not itself a failure).\n');
            end

            % Same seed, rebuilt from scratch: must reproduce bit-identically.
            objSeedA2 = InstanceSpace(classCaseDir, seedOptsA).build('stages', {'prelim', 'sifted', 'pilot'});
            testCase.verifyTrue(isequal(objSeedA.model.pilot.A, objSeedA2.model.pilot.A) && ...
                   isequal(objSeedA.model.sifted.selvars, objSeedA2.model.sifted.selvars), ...
                'Rebuilding with the same opts.general.seed changed PILOT.A/SIFTED.selvars -- reproducibility regressed.');
        end

        function testCloisterNonMeanCentredWarning(testCase)
            % #44: CLOISTER's correlation-contradiction filter assumes
            % mean-centred data; opts.auto.preproc=false or
            % opts.norm.flag=false skips PRELIM's Box-Cox+Z-score step
            % entirely, silently making the filter degenerate for any
            % naturally all-positive feature. runCloister should warn in
            % that case, and must not warn on the default (normalised) path.
            classCaseDir = testCase.CaseDir;
            baseOpts = testCase.BaseOpts;

            noNormOpts = baseOpts;
            noNormOpts.norm.flag = false;
            testCase.verifyWarning(@() InstanceSpace(classCaseDir, noNormOpts).build('stages', {'prelim', 'sifted', 'pilot', 'cloister'}), ...
                'ISA:InstanceSpace:cloisterNotMeanCentred', ...
                'opts.norm.flag=false should raise ISA:InstanceSpace:cloisterNotMeanCentred before CLOISTER runs.');

            lastwarn('');
            InstanceSpace(classCaseDir, baseOpts).build('stages', {'prelim', 'sifted', 'pilot', 'cloister'});
            [~, warnId2] = lastwarn();
            testCase.verifyNotEqual(warnId2, 'ISA:InstanceSpace:cloisterNotMeanCentred', ...
                'The default (normalised) path should not raise the cloisterNotMeanCentred warning.');
        end

        function testRequiredFieldsValidation(testCase)
            % #28: checkPrereq only confirms the prerequisite STAGE
            % completed; checkRequiredFields additionally confirms the
            % specific obj.model fields that stage needs are actually
            % present -- catches a hand-edited/incomplete model.mat that
            % completedStages alone would wave through.
            classCaseDir = testCase.CaseDir;
            baseOpts = testCase.BaseOpts;

            reqObj = InstanceSpace(classCaseDir, baseOpts).build('stages', {'prelim'});
            reqObj.model.data = rmfield(reqObj.model.data, 'Ybin'); % simulate an incomplete/corrupted model
            testCase.verifyError(@() reqObj.build('stages', {'sifted'}), 'ISA:InstanceSpace:missingField', ...
                ['build(''stages'',{''sifted''}) with obj.model.data.Ybin removed should raise ' ...
                 'ISA:InstanceSpace:missingField, not an opaque crash inside SIFTED.']);

            % Success path must be unaffected by the added check.
            reqObj2 = InstanceSpace(classCaseDir, baseOpts).build('stages', {'prelim', 'sifted'});
            testCase.verifyEqual(reqObj2.completedStages, {'prelim', 'sifted'}, ...
                'checkRequiredFields should not interfere with a normal successful build.');
        end

        function testTraceGoodBestMismatchGuard(testCase)
            % Folded-in finding from #28: TRACE's eval mode used
            % ngood/nbest as the real-vs-placeholder boundary for
            % out.good/out.best without ever checking they agree, relying
            % entirely on an unasserted cross-file convention.
            baseModel = testCase.BaseModel;
            testCase.assumeTrue(numel(baseModel.trace.good) >= 2, ...
                'This check needs at least 2 trained algorithms to desync good/best counts.');
            mismatchedTrace = baseModel.trace;
            mismatchedTrace.good = mismatchedTrace.good(1:end-1);
            testCase.verifyError(@() TRACE(baseModel.pilot.Z, baseModel.data.Ybin, baseModel.pythia.Yhat, ...
                baseModel.data.P, baseModel.data.beta, baseModel.data.algolabels, ...
                baseModel.opts.trace, mismatchedTrace), ...
                'ISA:TRACE:goodBestCountMismatch', ...
                'TRACE eval mode with mismatched trainedTrace.good/.best counts should raise ISA:TRACE:goodBestCountMismatch.');
        end

        function testPrelimTieBreakingConsistency(testCase)
            % Before #38, InstanceSpace.evaluateTestSet reimplemented
            % PRELIM's Ybin/Ybest/P/beta computation independently, using a
            % silent deterministic sort() with no tie-breaking logic at all
            % -- diverging from PRELIM's own random tie-break for instances
            % with more than one best-performing algorithm. #38 gave PRELIM
            % itself an eval mode (nargin==4) that shares this exact code
            % with training mode, so the two can no longer drift apart.
            % Verified directly here (not through a full build+explore
            % round trip) with a synthetic tie, seeding rng identically
            % before each call: PRELIM's tie-break is the only randomness
            % before that point in the function, so an identical seed must
            % give an identical pick if -- and only if -- both modes are
            % genuinely running the same code.
            Xsynth = [1 2; 3 4; 5 6; 7 8];
            Ysynth = [1 1 5; 2 3 4; 1 5 2; 6 2 1]; % row 1: algorithms 1 and 2 tied for best (both = 1)
            prelimTestOpts = struct('MaxPerf', false, 'AbsPerf', false, 'epsilon', 0.05, ...
                'betaThreshold', 0.55, 'auto', true, 'bound', true, 'norm', true, 'iqrMultiplier', 5);

            rng(7, 'twister');
            [~, ~, trainOut] = PRELIM(Xsynth, Ysynth, prelimTestOpts);
            testCase.verifyTrue(ismember(trainOut.P(1), [1, 2]), ...
                'Training-mode tie-break for the tied instance should pick one of the tied algorithms (1 or 2).');

            rng(7, 'twister');
            [~, ~, evalOut] = PRELIM(Xsynth, Ysynth, prelimTestOpts, trainOut);
            testCase.verifyEqual(evalOut.P(1), trainOut.P(1), ...
                ['Eval-mode PRELIM should break the same tie the same way training mode does when seeded ' ...
                 'identically -- #37''s drift meant these used to disagree.']);
        end
    end
end
