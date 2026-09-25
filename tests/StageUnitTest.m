classdef StageUnitTest < matlab.unittest.TestCase
% StageUnitTest  Direct tests of the standalone pipeline functions in core/
% on small synthetic data: the branches and input guards that a full
% pipeline run on the reference dataset does not reach (performance
% modes, degenerate inputs, precomputed solutions, error identifiers).
% Each test runs in seconds.

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
        Defaults    % ISAdefaults(struct()), for complete per-stage opts
    end

    properties (TestParameter)
        MaxPerf = {true, false};
        AbsPerf = {true, false};
        FilterType = {'Ftr', 'Ftr&AP', 'Ftr&Good', 'Ftr&AP&Good'};
        Classifier = {'tree', 'nb', 'linear'};
    end

    methods (TestClassSetup)
        function setup(testCase)
            root = testRepoRoot();
            addpath([root 'core'], [root 'utils'], [root 'output']);
            testCase.Defaults = ISAdefaults(struct());
            testCase.Defaults.pilot.verbose = false;
            testCase.Defaults.pythia.verbose = false;
        end
    end

    methods (TestMethodSetup)
        function closeFigures(testCase)
            testCase.addTeardown(@() close('all'));
        end
    end

    methods (Test)
        % ----------------------------------------------------------- PRELIM
        function testPrelimPerformanceModes(testCase, MaxPerf, AbsPerf)
            [X, Y] = syntheticMetadata(80, 4, 3);
            opts = prelimOpts(MaxPerf, AbsPerf);
            if AbsPerf
                opts.epsilon = median(Y(:)); % roughly half of the entries are good
            end
            [Xn, Yn, out] = PRELIM(X, Y, opts);
            testCase.verifySize(Xn, size(X));
            testCase.verifySize(Yn, size(Y));
            testCase.verifyClass(out.Ybin, 'logical');
            testCase.verifyTrue(all(out.P >= 1 & out.P <= size(Y, 2)));
            testCase.verifyEqual(out.numGoodAlgos, sum(out.Ybin, 2));
            if MaxPerf
                testCase.verifyEqual(out.Ybest, max(Y, [], 2));
            else
                testCase.verifyEqual(out.Ybest, min(Y, [], 2));
            end
            if ~AbsPerf
                testCase.verifyTrue(all(any(out.Ybin, 2)), ...
                    'With relative performance the best algorithm is always good.');
            end
            % The normalised columns have zero mean and unit variance.
            testCase.verifyEqual(mean(Xn), zeros(1, size(X, 2)), 'AbsTol', 1e-8);
            testCase.verifyEqual(std(Xn), ones(1, size(X, 2)), 'AbsTol', 1e-8);
        end

        function testPrelimWithoutPreprocessing(testCase)
            [X, Y] = syntheticMetadata(40, 3, 2);
            opts = prelimOpts(false, false);
            opts.auto = false;
            [Xn, Yn, out] = PRELIM(X, Y, opts);
            testCase.verifyEqual(Xn, X, 'opts.auto=false must leave X unchanged.');
            testCase.verifySize(Yn, size(Y));
            testCase.verifyEqual(out.lambdaX, zeros(1, 3));
        end

        function testPrelimGuards(testCase)
            [X, Y] = syntheticMetadata(10, 3, 2);
            testCase.verifyError(@() PRELIM(X(1:5,:), Y, prelimOpts(false, false)), ...
                'ISA:PRELIM:sizeMismatch');
            testCase.verifyError(@() PRELIM(X, Y, 'opts'), 'ISA:PRELIM:badOpts');
        end

        function testPrelimWarnsOnZeroBest(testCase, MaxPerf)
            % Relative performance divides by the best value; many exact
            % zeros there make the relative matrix meaningless.
            [X, Y] = syntheticMetadata(40, 3, 2);
            Y(1:10, :) = 0;
            testCase.verifyWarning(@() PRELIM(X, Y, prelimOpts(MaxPerf, false)), ...
                'ISA:PRELIM:manyZeroBest');
        end

        function testPrelimEvalModeClipsAndWarns(testCase)
            [X, Y] = syntheticMetadata(60, 3, 2);
            opts = prelimOpts(false, false);
            [~, ~, trained] = PRELIM(X, Y, opts);
            Xfar = X * 100; % far outside the training bounds
            [Xe, Ye, evalOut] = testCase.verifyWarning(@() PRELIM(Xfar, Y, opts, trained), ...
                'ISA:InstanceSpace:outOfDistribution');
            testCase.verifySize(Xe, size(X));
            testCase.verifyTrue(all(isfinite(Xe(:))) && all(isfinite(Ye(:))));
            testCase.verifyFalse(isfield(evalOut, 'lambdaX'), ...
                'Evaluation mode applies the trained parameters; it does not refit them.');
        end

        % ----------------------------------------------------------- FILTER
        function testFilterTypes(testCase, FilterType)
            [X, Y] = syntheticMetadata(50, 3, 2);
            X(26:50, :) = X(1:25, :) + 1e-6;       % near-duplicate of every instance
            Y(26:50, :) = Y(1:25, :);
            Ybin = true(size(Y));                   % every algorithm good everywhere
            opts = struct('mindistance', 0.01, 'type', FilterType);
            [redundant, isDissimilar, isVISA, unif] = FILTER(X, Y, Ybin, opts);
            testCase.verifyEqual(sum(redundant), 25, ...
                'Each near-duplicate should be removed under every condition here.');
            testCase.verifyFalse(any(isDissimilar(26:50)));
            testCase.verifyFalse(any(isVISA));
            testCase.verifyTrue(unif <= 1);
        end

        function testFilterKeepsCloseButDifferentInstances(testCase)
            [X, Y] = syntheticMetadata(20, 3, 2);
            X(11:20, :) = X(1:10, :);
            Y(11:20, :) = Y(1:10, :) + 10;         % very different performance
            Ybin = false(size(Y));                  % no instance is good everywhere
            for t = {'Ftr&AP', 'Ftr&Good', 'Ftr&AP&Good'}
                [redundant, ~, isVISA] = FILTER(X, Y, Ybin, struct('mindistance', 0.01, 'type', t{1}));
                testCase.verifyFalse(any(redundant), ['Nothing should be removed under ' t{1}]);
                testCase.verifyTrue(all(isVISA(11:20)), ['Close-but-kept instances are VISA under ' t{1}]);
            end
        end

        function testFilterGuards(testCase)
            [X, Y] = syntheticMetadata(10, 3, 2);
            testCase.verifyError(@() FILTER(X, Y, Y > 0, struct('mindistance', 0.1, 'type', 'X')), ...
                'ISA:FILTER:invalidType');
            same = repmat(X(1,:), 10, 1);           % every instance identical
            testCase.verifyWarning(@() FILTER(same, Y, Y > 0, struct('mindistance', 0.1, 'type', 'Ftr')), ...
                'ISA:FILTER:degenerateUniformity');
        end

        % ----------------------------------------------------------- SIFTED
        function testSiftedGuards(testCase)
            [X, Y] = syntheticMetadata(30, 4, 2);
            Ybin = Y > median(Y(:));
            labels = {'a', 'b', 'c', 'd'};
            opts = testCase.Defaults.sifted;
            testCase.verifyError(@() SIFTED(X(1:20,:), Y, Ybin, labels, opts), 'ISA:SIFTED:sizeMismatch');
            testCase.verifyError(@() SIFTED(X, Y, Ybin, labels(1:3), opts), 'ISA:SIFTED:labelMismatch');
            bad = opts; bad.dims = 4;
            testCase.verifyError(@() SIFTED(X, Y, Ybin, labels, bad), 'ISA:SIFTED:invalidDims');
            testCase.verifyError(@() SIFTED(X(:,1), Y, Ybin, labels(1), opts), 'ISA:SIFTED:tooFewFeatures');
        end

        function testSiftedKeepsFewFeatures(testCase)
            % Three or fewer features skip selection entirely.
            [X, Y] = syntheticMetadata(30, 3, 2);
            [Xs, out] = SIFTED(X, Y, Y > median(Y(:)), {'a', 'b', 'c'}, testCase.Defaults.sifted);
            testCase.verifyEqual(out.selvars, 1:3);
            testCase.verifyEqual(Xs, X);
        end

        % ------------------------------------------------------------ PILOT
        function testPilotGuards(testCase)
            [X, Y] = syntheticMetadata(30, 4, 2);
            labels = {'a', 'b', 'c', 'd'};
            opts = testCase.Defaults.pilot;
            bad = opts; bad.method = 'pca';
            testCase.verifyError(@() PILOT(X, Y, labels, bad), 'ISA:PILOT:invalidMethod');
            bad = opts; bad.dims = 5;
            testCase.verifyError(@() PILOT(X, Y, labels, bad), 'ISA:PILOT:invalidDims');
            bad = opts; bad.alpha = -1;
            testCase.verifyError(@() PILOT(X, Y, labels, bad), 'ISA:PILOT:invalidAlpha');
        end

        function testPilotRankDeficientFallsBack(testCase)
            [X, Y] = syntheticMetadata(30, 4, 2);
            X(:, 4) = X(:, 1) + X(:, 2);            % rank deficient
            opts = testCase.Defaults.pilot;
            opts.analytic = true;
            opts.ntries = 2;
            out = testCase.verifyWarning(@() PILOT(X, Y, {'a', 'b', 'c', 'd'}, opts), ...
                'ISA:PILOT:rankDeficient');
            testCase.verifySize(out.Z, [30 2]);
            testCase.verifyTrue(isfield(out, 'alpha'), 'The fallback is the numerical solver.');
        end

        function testPilotReusesSolutions(testCase)
            [X, Y] = syntheticMetadata(30, 3, 2);
            labels = {'a', 'b', 'c'};
            opts = testCase.Defaults.pilot;
            opts.ntries = 2;
            first = PILOT(X, Y, labels, opts);

            % Given starting points reproduce the same trials.
            withX0 = opts; withX0.X0 = first.X0;
            again = PILOT(X, Y, labels, withX0);
            testCase.verifyEqual(again.Z, first.Z, 'AbsTol', 1e-10);

            % A precomputed solution skips the optimisation.
            [~, best] = max(first.perf);
            precalc = opts; precalc.precalcAlpha = first.alpha(:, best);
            reused = PILOT(X, Y, labels, precalc);
            testCase.verifyEqual(reused.Z, first.Z, 'AbsTol', 1e-10);
        end

        function testPilotViewpointGuardsAndStartPoints(testCase)
            [X, Y] = syntheticMetadata(30, 3, 2);
            opts = testCase.Defaults.pilot;
            opts.ntries = 2;
            testCase.verifyError(@() PILOTviewpoint(X(:,1:2), Y, opts), 'ISA:PILOTviewpoint:not3D');
            testCase.verifyError(@() PILOTviewpoint(X(1:10,:), Y, opts), 'ISA:PILOTviewpoint:sizeMismatch');
            opts.X0 = 2*rand(2*3 + 2*size(Y, 2), 2) - 1;
            vp = PILOTviewpoint(X, Y, opts);
            testCase.verifySize(vp.A{1}, [2 3]);
            testCase.verifyTrue(isfinite(vp.azimuth(1)) && isfinite(vp.elevation(1)));
        end

        % --------------------------------------------------------- CLOISTER
        function testCloisterFeatureLimit(testCase)
            [X, ~] = syntheticMetadata(40, 4, 2);
            opts = testCase.Defaults.cloister;
            opts.maxFeatures = 3;
            A2 = randn(2, 4);
            out = testCase.verifyWarning(@() CLOISTER(X, A2, opts), 'ISA:CLOISTER:tooManyFeatures');
            testCase.verifyEqual(out.Zedge(1,:), out.Zedge(end,:), 'The 2D fallback is a closed polygon.');
            testCase.verifyEqual(out.Zecorr, out.Zedge);

            A3 = randn(3, 4);
            out3 = testCase.verifyWarning(@() CLOISTER(X, A3, opts), 'ISA:CLOISTER:tooManyFeatures');
            testCase.verifyEqual(size(out3.ZedgeFaces, 2), 3, 'The 3D fallback is a triangulated hull.');
        end

        function testCloisterCoplanar3D(testCase)
            % A 3D projection of two features puts every corner on one
            % plane: there is no volumetric hull, so CLOISTER must return
            % the flat polygon, triangulated, instead of failing in
            % convhull.
            [X, ~] = syntheticMetadata(40, 2, 2);
            out = CLOISTER(X, randn(3, 2), testCase.Defaults.cloister);
            testCase.verifyEqual(size(out.Zedge, 2), 3);
            testCase.verifyEqual(size(out.ZedgeFaces, 2), 3);
            testCase.verifyEqual(size(out.ZedgeFaces, 1), size(out.Zedge, 1) - 2, ...
                'A fan triangulation of an n-gon has n-2 triangles.');
            testCase.verifyLessThanOrEqual(max(out.ZedgeFaces(:)), size(out.Zedge, 1));
        end

        function testCloisterStrictThresholdFallsBack(testCase)
            % Perfectly correlated features and a zero threshold discard all
            % but two corners, too few for a hull: Zecorr falls back to Zedge.
            base = linspace(-1, 1, 40)';
            X = [base, 2*base, 3*base] + 1e-3*randn(40, 3);
            opts = testCase.Defaults.cloister;
            opts.corrThreshold = 0;
            out = CLOISTER(X, randn(2, 3), opts);
            testCase.verifyEqual(out.Zecorr, out.Zedge);
        end

        % ------------------------------------------------------------ TRACE
        function testTraceWithoutPredictions(testCase)
            [Z, Ybin, P, beta] = syntheticSpace();
            opts = testCase.Defaults.trace;
            out = testCase.verifyWarning(@() TRACE(Z, Ybin, [], P, beta, {'a', 'b'}, opts), ...
                'ISA:TRACE3:noPYTHIA');
            testCase.verifyNumElements(out.good, 2);
            testCase.verifyGreaterThan(out.good{1}.measure, 0, ...
                'A dense cluster of good instances should give a footprint.');
        end

        function testTraceLegacyFallsBackIn3D(testCase)
            [Z, Ybin, P, beta] = syntheticSpace();
            Z = [Z, randn(size(Z, 1), 1)];
            opts = testCase.Defaults.trace;
            opts.method = 'legacy';
            out = testCase.verifyWarning(@() TRACE(Z, Ybin, Ybin, P, beta, {'a', 'b'}, opts), ...
                'ISA:TRACE:legacyNo3D');
            testCase.verifyEqual(out.space.measureLabel, 'Volume');
        end

        function testTraceLegacy2DEvaluation(testCase)
            [Z, Ybin, P, beta] = syntheticSpace();
            opts = testCase.Defaults.trace;
            opts.method = 'legacy';
            opts.contra = true;
            trained = TRACE(Z, Ybin, Ybin, P, beta, {'a', 'b'}, opts);
            testCase.verifyNumElements(trained.best, 2);
            evaluated = TRACE(Z, Ybin, Ybin, P, beta, {'a', 'b'}, opts, trained);
            testCase.verifySize(evaluated.summary, size(trained.summary));
        end

        % ----------------------------------------------------------- PYTHIA
        function testPythiaSingleParameterClassifiers(testCase, Classifier)
            [Z, Ybin, ~, ~, Y] = syntheticSpace();
            opts = testCase.Defaults.pythia;
            opts.classifier = Classifier;
            opts.nTuningIter = 3;
            opts.kFold = 3;
            out = PYTHIA(Z, Y, Ybin, min(Y, [], 2), {'a', 'b'}, opts);
            testCase.verifyNumElements(out.classifiers, 2);
            testCase.verifySize(out.summary, [2+3, 10], ...
                'One-parameter classifiers have 10 summary columns.');
            evaluated = PYTHIA(Z, Y, Ybin, min(Y, [], 2), {'a', 'b'}, opts, out);
            testCase.verifySize(evaluated.Yhat, size(Ybin));
        end

        function testPythiaDegenerateLabels(testCase)
            % An algorithm that is good everywhere gets a constant
            % predictor instead of a classifier, in training and evaluation.
            [Z, Ybin, ~, ~, Y] = syntheticSpace();
            Ybin(:, 2) = true;
            opts = testCase.Defaults.pythia;
            opts.verbose = true;
            opts.nTuningIter = 2;
            opts.kFold = 3;
            out = testCase.verifyWarning(@() PYTHIA(Z, Y, Ybin, min(Y, [], 2), {'a', 'b'}, opts), ...
                'ISA:PYTHIA:degenerateLabel');
            testCase.verifyTrue(all(out.Yhat(:, 2)));
            evaluated = PYTHIA(Z, Y, Ybin, min(Y, [], 2), {'a', 'b'}, opts, out);
            testCase.verifyTrue(all(evaluated.Yhat(:, 2)));
        end

        function testPythiaConstantWeights(testCase)
            [Z, Ybin] = syntheticSpace();
            Y = ones(size(Ybin));                   % constant performance: no usable weights
            opts = testCase.Defaults.pythia;
            opts.useweights = true;
            opts.nTuningIter = 2;
            opts.kFold = 3;
            testCase.verifyWarning(@() PYTHIA(Z, Y, Ybin, min(Y, [], 2), {'a', 'b'}, opts), ...
                'ISA:PYTHIA:degenerateWeights');
        end

        function testPythiaGuards(testCase)
            [Z, Ybin, ~, ~, Y] = syntheticSpace();
            Ybest = min(Y, [], 2);
            opts = testCase.Defaults.pythia;
            bad = opts; bad.nTuningIter = 0;
            testCase.verifyError(@() PYTHIA(Z, Y, Ybin, Ybest, {'a', 'b'}, bad), ...
                'ISA:PYTHIA:invalidNTuningIter');
            bad = opts; bad.tuning = 'none'; bad.params = [];
            testCase.verifyError(@() PYTHIA(Z, Y, Ybin, Ybest, {'a', 'b'}, bad), ...
                'ISA:PYTHIA:noParamsForNoneTuning');
            bad = opts; bad.uselibsvm = true; bad.skip = true;
            testCase.verifyWarning(@() PYTHIA(Z, Y, Ybin, Ybest, {'a', 'b'}, bad), ...
                'ISA:PYTHIA:libsvmDeprecated');
            testCase.verifyError(@() PYTHIA(Z, Y, Ybin, Ybest, {'a', 'b'}, opts, struct('mu', 0)), ...
                'ISA:PYTHIA:noClassifier');
        end

        % ---------------------------------------------------- ISArecallView
        function testRecallViewGuards(testCase)
            testCase.verifyError(@() ISArecallView(5), 'ISA:ISArecallView:notAFigure');
            fig = figure('Visible', 'off');
            testCase.verifyError(@() ISArecallView(fig), 'ISA:ISArecallView:noStoredViewpoint');
        end

        function testRecallViewRestoresAngle(testCase)
            fig = figure('Visible', 'off');
            vp = struct('groups', {{[1 2], 3}}, 'azimuth', [0.5; 1.0], 'elevation', [0.2; 0.4]);
            fig.UserData = struct('isaViewpoint', vp);
            testCase.verifyError(@() ISArecallView(fig), 'ISA:ISArecallView:noAxes');

            ax = axes(fig);
            scatter3(ax, rand(10,1), rand(10,1), rand(10,1));
            ISArecallView(fig, 3);                  % algorithm 3 is in group 2
            [az, el] = view(ax);
            testCase.verifyEqual([az el], rad2deg([1.0 0.4]), 'AbsTol', 1e-6);
            ISArecallView(fig);                     % global viewpoint
            [az, el] = view(ax);
            testCase.verifyEqual([az el], rad2deg([0.5 0.2]), 'AbsTol', 1e-6);
        end
    end
end

% =========================================================================
function [X, Y] = syntheticMetadata(ninst, nfeats, nalgos)
% Positive features and performance with a linear dependence between them.
rng(11, 'twister');
X = rand(ninst, nfeats) + 0.1;
W = rand(nfeats, nalgos);
Y = X*W + 0.1*rand(ninst, nalgos);
end

function opts = prelimOpts(maxPerf, absPerf)
opts = struct('MaxPerf', maxPerf, 'AbsPerf', absPerf, 'epsilon', 0.1, ...
              'betaThreshold', 0.55, 'auto', true, 'bound', true, ...
              'norm', true, 'iqrMultiplier', 5);
end

function [Z, Ybin, P, beta, Y] = syntheticSpace()
% Two separated clusters in 2D: algorithm 1 is good in the first, and
% algorithm 2 in the second.
rng(5, 'twister');
n = 60;
Z = [randn(n, 2)*0.5; randn(n, 2)*0.5 + 4];
Ybin = [[true(n, 1); false(n, 1)], [false(n, 1); true(n, 1)]];
Y = double(~Ybin) + 0.1*rand(2*n, 2);       % a cost: good means low
[~, P] = min(Y, [], 2);
beta = any(Ybin, 2);
end
