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
            % Absolute, not './test/data/': matlab.unittest does not
            % guarantee the working directory during test execution is
            % the directory test_integration.m was launched from (see
            % testRepoRoot.m).
            rootdir = [testRepoRoot() 'test/data/'];
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

            % PR #51 review: runSifted() unconditionally read-then-writes
            % obj.model.featsel.idx (obj.model.featsel.idx =
            % obj.model.featsel.idx(obj.model.sifted.selvars)), a genuine
            % prerequisite the 'sifted' contract had omitted -- a model
            % missing it crashed opaquely inside runSifted instead of via
            % the intended clear error.
            reqObj3 = InstanceSpace(classCaseDir, baseOpts).build('stages', {'prelim'});
            reqObj3.model = rmfield(reqObj3.model, 'featsel');
            testCase.verifyError(@() reqObj3.build('stages', {'sifted'}), 'ISA:InstanceSpace:missingField', ...
                ['build(''stages'',{''sifted''}) with obj.model.featsel.idx removed should raise ' ...
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

        function testBoundaryDisplay2D(testCase)
            % #32 (option (a), 2D only): CLOISTER's boundary was computed
            % but never rendered in automated output or exposed via
            % plot(). BaseModel is a 2D build (testDefaultOpts default),
            % so both the scriptpng.m PNG and plot('boundary') should now
            % work without error.
            baseModel = testCase.BaseModel;
            testCase.assumeTrue(isfield(baseModel, 'cloist'), ...
                'This check needs a model with CLOISTER built.');
            pngPath = [testCase.CaseDir 'distribution_boundary.png'];
            testCase.verifyTrue(isfile(pngPath), ...
                'scriptpng.m should have written distribution_boundary.png for a 2D build with CLOISTER run.');

            obj = InstanceSpace(testCase.CaseDir, testCase.BaseOpts);
            obj.model = baseModel;
            fig = figure('Visible', 'off');
            testCase.addTeardown(@() close(fig));
            % No try/catch, no verifyWarningFree: an uncaught error here
            % fails this test automatically, and headless CI can
            % legitimately warn (e.g. graphics-acceleration-unavailable)
            % without that being a real failure.
            obj.plot('boundary');
        end

        function testBoundaryDisplayRejects3D(testCase)
            % The 'boundary' view must refuse a 3D projection rather than
            % draw an inaccurate boundary: CLOISTER's Zedge/Zecorr are
            % computed via a 2D-only convex hull (core/CLOISTER.m) even
            % when the projection itself is 3D (#32's own scope note).
            % Synthesised directly (no full 3D build) since only
            % obj.model.pilot.Z's column count and obj.model.cloist's
            % presence are checked before this error is raised.
            obj = InstanceSpace(testCase.CaseDir, testCase.BaseOpts);
            obj.model.pilot.Z = zeros(5, 3);
            obj.model.cloist = struct('Zedge', zeros(5, 2));
            testCase.verifyError(@() obj.plot('boundary'), 'ISA:InstanceSpace:boundaryNot3D', ...
                'plot(''boundary'') on a 3D model should raise ISA:InstanceSpace:boundaryNot3D, not draw an inaccurate boundary.');
        end

        function testTraceAlphaBoundaryMultiRegion(testCase)
            % #31: traceAlphaBoundary previously built one adjacency graph
            % from boundaryFacets(poly)'s combined (all-regions) edge list
            % and stopped as soon as it ran out of same-loop neighbours,
            % silently returning only the first region traced for a
            % disconnected/multi-region alpha shape. Reproduced directly
            % with two well-separated point clusters (#31's own
            % acceptance criteria step 1: reproduce the multi-region
            % case) and confirmed here that every region's vertices come
            % back now, not just the first.
            scriptfcn; % injects traceAlphaBoundary into this function's workspace
            cluster1 = [0 0; 1 0; 0 1; 1 1];
            cluster2 = [20 20; 21 20; 20 21; 21 21];
            pts = [cluster1; cluster2];
            as = alphaShape(pts, 2); % alpha << inter-cluster gap, keeps the two regions disconnected
            testCase.assumeTrue(numRegions(as) >= 2, ...
                'This check needs a genuinely multi-region alphaShape to reproduce #31.');

            verts = traceAlphaBoundary(as); %#ok<NODEF> -- injected by scriptfcn above
            testCase.verifyTrue(any(all(isnan(verts), 2)), ...
                'traceAlphaBoundary should separate multiple regions with a NaN row, not silently drop all but one.');
            testCase.verifyTrue(any(verts(:,1) < 5) && any(verts(:,1) > 15), ...
                'traceAlphaBoundary should include vertices from BOTH well-separated clusters, not just the first region traced.');
        end

        function testPrelimEvalModeAlgoAlignmentAfterPruning(testCase)
            % PR #51 review: InstanceSpace.runPrelim removes any algorithm
            % with no good instances from data.Yraw/Y/Ybin/algolabels, but
            % left trainedPrelim.lambdaY/muY/sigmaY (per-algorithm,
            % fit before that pruning) untouched. PRELIM's eval mode
            % derives modelalgos from numel(trainedPrelim.lambdaY), so an
            % unpruned lambdaY both over-counts modelalgos against the
            % actually-reconciled eval-time Y (indexing past it, or
            % misapplying a pruned algorithm's transform to an unrelated
            % new algorithm's column) and, if the pruned algorithm wasn't
            % last, misaligns every surviving lambda/mu/sigma positioned
            % after it. Reproduced directly against PRELIM rather than
            % the full pipeline: forcing InstanceSpace.build() to prune a
            % real algorithm needs data engineered so one is never
            % "good," which is fragile against the bundled reference
            % dataset; simulating the post-prune trainedPrelim state
            % directly is equivalent and deterministic. Middle algorithm
            % (not the last) pruned deliberately, since that's what
            % exposes the positional-misalignment half of the bug, not
            % just the over-counting half.
            opts = struct('MaxPerf', false, 'AbsPerf', true, 'epsilon', 10, ...
                'betaThreshold', 0.55, 'auto', true, 'bound', true, ...
                'norm', true, 'iqrMultiplier', 5);
            rng(1, 'twister');
            Xtrain = rand(40, 3) + 1;
            Ytrain = rand(40, 3) + 1;
            [~, ~, trainedPrelim] = PRELIM(Xtrain, Ytrain, opts);

            keep = [true false true]; % algorithm 2 of 3 "pruned"
            trainedPrelim.lambdaY = trainedPrelim.lambdaY(keep);
            trainedPrelim.muY     = trainedPrelim.muY(keep);
            trainedPrelim.sigmaY  = trainedPrelim.sigmaY(keep);

            Xtest = rand(10, 3) + 1;
            Ytest = rand(10, 3) + 1;
            Ytest = Ytest(:, keep); % as INIT.m/evaluateTestSet would reconcile it

            [~, YOut, ~] = PRELIM(Xtest, Ytest, opts, trainedPrelim);
            testCase.verifyEqual(size(YOut, 2), 2, ...
                'PRELIM eval mode should not resize/error on Y when trainedPrelim''s per-algorithm fields are pruned to match the reconciled algorithm count.');
            testCase.verifyTrue(all(isfinite(YOut(:))), ...
                'PRELIM eval-mode normalisation produced non-finite values after algorithm pruning -- likely misapplied a wrong algorithm''s Box-Cox/Z-score transform.');
        end
    end
end
