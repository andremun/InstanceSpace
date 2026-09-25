classdef InstanceSpaceTest < matlab.unittest.TestCase
% InstanceSpaceTest  The InstanceSpace class's error paths and plot views,
% instance sources, the web output, file-index subsetting, and explore()
% on test metadata whose algorithms or features differ from training.
% One model is built once, from the reference data with a 'source' column
% added, and shared by every test.

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
        Opts
        Built       % fully built InstanceSpace object
    end

    methods (TestClassSetup)
        function buildModel(testCase)
            root = [testRepoRoot() 'test/data/'];
            testCase.CaseDir = [root 'instance_space_api/'];
            if ~isfolder(testCase.CaseDir), mkdir(testCase.CaseDir); end
            T = readtable([root 'metadata.csv']);
            suites = {'suiteA'; 'suiteB'};
            T.source = suites(mod((1:height(T))', 2) + 1);
            writetable(T, [testCase.CaseDir 'metadata.csv']);
            copyfile([root 'metadata_test.csv'], [testCase.CaseDir 'metadata_test.csv']);

            testCase.Opts = testDefaultOpts();
            testCase.Opts.outputs.web = true;     % covers scriptweb
            testCase.Built = InstanceSpace(testCase.CaseDir, testCase.Opts).build();
        end
    end

    methods (TestMethodSetup)
        function invisibleFigure(testCase)
            fig = figure('Visible', 'off');
            testCase.addTeardown(@() close(fig));
        end
    end

    methods (Test)
        function testSourcesAndWebOutputs(testCase)
            m = testCase.Built.model;
            testCase.verifyClass(m.data.S, 'categorical');
            testCase.verifyEqual(sort(categories(m.data.S)), {'suiteA'; 'suiteB'});
            testCase.verifyTrue(isfile([testCase.CaseDir 'distribution_sources.png']), ...
                'scriptpng should draw the instance sources.');
            for f = {'color_table.csv', 'feature_raw_color.csv', 'algorithm_process_single_color.csv', 'good_algos_color.csv'}
                testCase.verifyTrue(isfile([testCase.CaseDir f{1}]), ...
                    ['opts.outputs.web=true should write ' f{1} '.']);
            end
            testCase.Built.plot('sources');
        end

        function testPlotViews(testCase)
            obj = testCase.Built;
            obj.plot('portfolio');
            obj.plot('good', 1);
            obj.plot('footprint', 1);
            obj.plot('boundary');
            nalgos = numel(obj.model.data.algolabels);
            testCase.verifyError(@() obj.plot('coverage'), 'ISA:InstanceSpace:unknownView');
            testCase.verifyError(@() obj.plot('good'), 'ISA:InstanceSpace:missingAlgoIdx');
            testCase.verifyError(@() obj.plot('footprint', nalgos + 1), 'ISA:InstanceSpace:badAlgoIdx');
            testCase.verifyError(@() obj.plot('good', 0), 'ISA:InstanceSpace:badAlgoIdx');

            noBound = obj;
            noBound.model = rmfield(noBound.model, 'cloist');
            testCase.verifyError(@() noBound.plot('boundary'), 'ISA:InstanceSpace:noCloister');
            noSource = obj;
            noSource.model.data = rmfield(noSource.model.data, 'S');
            testCase.verifyError(@() noSource.plot('sources'), 'ISA:InstanceSpace:noSources');
        end

        function testUnbuiltObjectErrors(testCase)
            fresh = InstanceSpace(testCase.CaseDir, testCase.Opts);
            testCase.verifyError(@() fresh.plot('portfolio'), 'ISA:InstanceSpace:notBuilt');
            testCase.verifyError(@() fresh.save(), 'ISA:InstanceSpace:notBuilt');
            testCase.verifyError(@() fresh.explore(testCase.CaseDir), 'ISA:InstanceSpace:notBuilt');
            testCase.verifyEqual(fresh.getResults(), struct(), 'An unbuilt object has an empty model.');
            testCase.verifyError(@() testCase.Built.getResults(1), 'ISA:InstanceSpace:badResultIndex');
        end

        function testMissingFiles(testCase)
            empty = [testCase.CaseDir 'empty'];     % no trailing separator on purpose
            if ~isfolder(empty), mkdir(empty); end
            testCase.verifyError(@() InstanceSpace(empty), 'ISA:InstanceSpace:missingData');
            testCase.verifyError(@() InstanceSpace.load(empty), 'ISA:InstanceSpace:missingModel');
            testCase.verifyError(@() testCase.Built.explore(empty), 'ISA:InstanceSpace:missingTestData');
        end

        function testOptionsFromJsonFile(testCase)
            % Without an opts argument the constructor reads options.json.
            opts = struct('pilot', struct('dims', 3), 'perf', struct('epsilon', 0.3));
            fid = fopen([testCase.CaseDir 'options.json'], 'w');
            fprintf(fid, '%s', jsonencode(opts));
            fclose(fid);
            testCase.addTeardown(@() delete([testCase.CaseDir 'options.json']));
            obj = InstanceSpace(testCase.CaseDir);
            testCase.verifyEqual(obj.opts.pilot.dims, 3);
            testCase.verifyEqual(obj.opts.perf.epsilon, 0.3);
            testCase.verifyEqual(obj.opts.pilot.ntries, 10, 'Unset fields take their defaults.');
        end

        function testFileIndexSubset(testCase)
            idxFile = [testCase.CaseDir 'subset_index.csv'];
            writetable(table((1:2:100)', 'VariableNames', {'index'}), idxFile);
            opts = testCase.Opts;
            opts.selvars.fileidxflag = true;
            opts.selvars.fileidx = idxFile;
            obj = InstanceSpace(testCase.CaseDir, opts).build('stages', {'prelim'});
            testCase.verifyEqual(size(obj.model.data.X, 1), 50, ...
                'Only the instances listed in opts.selvars.fileidx should be kept.');
            testCase.verifyEqual(numel(obj.model.data.S), 50);
        end

        function testExploreReconcilesAlgorithms(testCase)
            % metadata_test.csv without one trained algorithm and with one
            % new algorithm. The missing algorithm must not be scored (#58);
            % the new one has no classifier; both are NaN in the summary.
            m = testCase.Built.model;
            testDir = [testCase.CaseDir 'explore_algos/'];
            if ~isfolder(testDir), mkdir(testDir); end
            T = readtable([testCase.CaseDir 'metadata_test.csv']);
            dropped = ['algo_' m.data.algolabels{1}];
            T.algo_NEW = T.(['algo_' m.data.algolabels{2}]);
            T.(dropped) = [];
            writetable(T, [testDir 'metadata_test.csv']);

            % Outputs off: this test is about the evaluation itself.
            obj = testCase.Built;
            obj.model.opts.outputs.csv = false;
            obj.model.opts.outputs.png = false;
            obj.model.opts.outputs.web = false;
            obj = obj.explore(testDir);
            res = obj.getResults(1);
            nalgos = numel(m.data.algolabels);
            ntest = numel(res.data.algolabels);
            testCase.verifyEqual(res.data.algolabels(1:nalgos), m.data.algolabels, ...
                'Trained algorithms keep their trained column positions.');
            testCase.verifyEqual(res.data.algolabels{end}, 'NEW', 'New algorithms are appended.');
            testCase.verifySize(res.pythia.Yhat, [height(T), ntest]);
            testCase.verifyTrue(all(isnan(res.data.Yraw(:, 1))), ...
                'A trained algorithm missing from the test data has no performance values.');
            testCase.verifyTrue(isnan(res.pythia.accuracy(1)), ...
                'A trained algorithm missing from the test data must not be scored (#58).');
            testCase.verifyTrue(isnan(res.pythia.accuracy(end)), ...
                'An algorithm with no trained classifier has no accuracy.');
            testCase.verifyTrue(all(isfinite(res.pythia.accuracy(2:nalgos))), ...
                'Algorithms present in both files are scored.');
            testCase.verifyNumElements(res.trace.good, ntest);
        end

        function testExploreRejectsFeatureMismatch(testCase)
            T = readtable([testCase.CaseDir 'metadata_test.csv']);
            names = T.Properties.VariableNames;
            feats = find(startsWith(names, 'feature_'));

            swapDir = [testCase.CaseDir 'explore_swapped/'];
            if ~isfolder(swapDir), mkdir(swapDir); end
            order = 1:numel(names);
            order(feats([1 2])) = feats([2 1]);
            writetable(T(:, order), [swapDir 'metadata_test.csv']);
            testCase.verifyError(@() testCase.Built.explore(swapDir), ...
                'ISA:InstanceSpace:featureOrderMismatch');

            dropDir = [testCase.CaseDir 'explore_dropped/'];
            if ~isfolder(dropDir), mkdir(dropDir); end
            T(:, feats(1)) = [];
            writetable(T, [dropDir 'metadata_test.csv']);
            testCase.verifyError(@() testCase.Built.explore(dropDir), ...
                'ISA:InstanceSpace:featureCountMismatch');
        end
    end
end
