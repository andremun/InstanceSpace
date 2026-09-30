classdef PortfolioPruningTest < matlab.unittest.TestCase
% PortfolioPruningTest  Keep preprocessing aligned with retained algorithms.
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
        MaxPerf = {false, true};
        AbsPerf = {false, true};
        Normalize = {false, true};
        DropColumn = {1, 2};
    end
    methods (Test)
        function testRetainedPortfolio(testCase, MaxPerf, AbsPerf, Normalize, DropColumn)
            root = testRepoRoot();
            addpath(root, [root 'core'], [root 'utils'], [root 'output']);
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            state = rng;
            testCase.addTeardown(@() rng(state));
            X = [(1:8)', [3; 1; 6; 2; 8; 4; 7; 5]];
            retained = [1 2; 2 1; 1 1; 3 8; 8 3; 7 9; 9 7; 2 3];
            dropped = 10*ones(8, 1);
            if AbsPerf
                % The removed algorithm wins this row, but none is good.
                dropped(6) = 6;
            end
            if MaxPerf
                retained = 20-retained;
                dropped = 20-dropped;
            end
            Y = [retained(:,1:DropColumn-1), dropped, retained(:,DropColumn:end)];
            labels = {'algo_a', 'algo_b', 'algo_c'};
            T = array2table([X, Y], 'VariableNames', [{'feature_x', 'feature_y'}, labels]);
            T.instances = cellstr("instance_" + (1:8)');
            writetable(T, fullfile(fixture.Folder, 'metadata.csv'));
            opts = ISAdefaults(testDefaultOpts());
            opts.general.verbose = false;
            opts.auto.preproc = Normalize;
            opts.perf.MaxPerf = MaxPerf;
            opts.perf.AbsPerf = AbsPerf;
            opts.perf.betaThreshold = 0.4;
            if AbsPerf
                opts.perf.epsilon = 5;
                if MaxPerf, opts.perf.epsilon = 15; end
            else
                opts.perf.epsilon = 0.05;
            end
            obj = InstanceSpace(fixture.Folder, opts).build('stages', {'prelim'});
            data = obj.model.data;
            p = opts.perf;
            p.auto = Normalize;
            p.bound = opts.bound.flag;
            p.norm = opts.norm.flag;
            p.iqrMultiplier = opts.prelim.iqrMultiplier;
            rng(opts.general.seed, 'twister');
            [expectedX, expectedY, expected] = PRELIM(X, retained, p);
            expectedLabels = {'a', 'b', 'c'};
            expectedLabels(DropColumn) = [];
            testCase.verifyEqual(data.algolabels, expectedLabels);
            testCase.verifyTrue(all(data.P >= 1 & data.P <= 2));
            testCase.verifyEqual(data.Yraw, retained);
            testCase.verifyEqual(data.X, expectedX, 'AbsTol', 1e-10);
            testCase.verifyEqual(data.Y, expectedY, 'AbsTol', 1e-10);
            fields = {'P', 'Ybest', 'Ybin', 'numGoodAlgos', 'beta'};
            for k = 1:numel(fields)
                field = fields{k};
                testCase.verifyEqual(data.(field), expected.(field), field);
                testCase.verifyEqual(obj.model.prelim.(field), expected.(field), field);
            end
            testCase.verifyEqual(obj.model.prelim.lambdaY, expected.lambdaY, 'AbsTol', 1e-10);
            testCase.verifyEqual(obj.model.prelim.muY, expected.muY, 'AbsTol', 1e-10);
            testCase.verifyEqual(obj.model.prelim.sigmaY, expected.sigmaY, 'AbsTol', 1e-10);
        end
    end
end
