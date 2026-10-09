classdef StateReviewTest < matlab.unittest.TestCase
% StateReviewTest  Regression cases for the September 2026 review.
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
        Folder
        Opts
    end
    methods (TestMethodSetup)
        function setup(tc)
            root = testRepoRoot(); addpath(root,[root 'core'],[root 'utils'],[root 'output']);
            fixture = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            tc.Folder = fixture.Folder;
            rng(13); X = rand(30,4); Y = rand(30,2);
            T = array2table([X Y], 'VariableNames', {'feature_a','feature_b','feature_c','feature_d','algo_a','algo_b'});
            T.instances = cellstr("i"+(1:30)');
            writetable(T, fullfile(tc.Folder,'metadata.csv'));
            tc.Opts = ISAdefaults(struct()); tc.Opts.sifted.flag = false;
            tc.Opts.pilot.analytic = true; tc.Opts.general.verbose = false;
            tc.Opts.outputs.csv = false; tc.Opts.outputs.png = false; tc.Opts.outputs.fig = false;
        end
    end
    methods (Test)
        function testPLSExploreMatchesTraining(tc)
            copyfile(fullfile(tc.Folder,'metadata.csv'),fullfile(tc.Folder,'metadata_test.csv'));
            for normalize = [false true]
                opts = tc.Opts; opts.auto.preproc = normalize; opts.pilot.method = 'pls';
                opts.selvars.smallscaleflag = true; opts.selvars.smallscale = .6; opts.pythia.skip = true;
                obj = InstanceSpace(tc.Folder,opts).build();
                obj = obj.explore(tc.Folder);
                result = obj.getResults(1);
                [found,rows] = ismember(obj.model.data.instlabels,result.data.instlabels);
                tc.verifyTrue(all(found));
                tc.verifyEqual(result.pilot.Z(rows,:),obj.model.pilot.Z,'AbsTol',1e-8);
            end
        end
        function testOwnedPoolClosesOnError(tc)
            tc.assumeTrue(isempty(gcp('nocreate')));
            opts = tc.Opts; opts.general.parallel = true; opts.general.ncores = 2;
            obj = InstanceSpace(tc.Folder,opts);
            tc.verifyError(@() obj.build('stages',{'prelim'},'onStage', ...
                @(varargin) error('ISA:test:stop','Stop after preprocessing.')),'ISA:test:stop');
            tc.verifyTrue(isempty(gcp('nocreate')));
        end
        function testStringPathAndSubsetValidation(tc)
            obj = InstanceSpace(string(tc.Folder),tc.Opts).build('stages',{'prelim'});
            tc.verifyEqual(obj.rootdir,[tc.Folder filesep]);
            opts = tc.Opts; opts.selvars.type = "ftr"; opts.selvars.densityflag = true;
            obj = InstanceSpace(string(tc.Folder),opts).build('stages',{'prelim'});
            tc.verifyTrue(obj.model.prelim.bydensity);
            opts.selvars.fileidxflag = true;
            tc.verifyError(@() InstanceSpace(tc.Folder,opts),'ISA:ISAvalidateOpts:subsetConflict');
            opts.selvars.densityflag = false; opts.selvars.fileidx = fullfile(tc.Folder,'absent.csv');
            tc.verifyError(@() InstanceSpace(tc.Folder,opts).build('stages',{'prelim'}),'ISA:InstanceSpace:missingIndexFile');
            opts.selvars.fileidx = fullfile(tc.Folder,'indices.csv');
            writetable(table([1;31]),opts.selvars.fileidx);
            tc.verifyError(@() InstanceSpace(tc.Folder,opts).build('stages',{'prelim'}),'ISA:InstanceSpace:badIndices');
            opts = ISAvalidateOpts(struct('pythia',struct('tuning','NONE')));
            tc.verifyEqual(opts.pythia.tuning,'none');
        end
        function testDroppedFeatureReplay(tc)
            file = fullfile(tc.Folder,'metadata.csv'); T = readtable(file);
            T.feature_b(:) = NaN; writetable(T,file);
            obj = InstanceSpace(tc.Folder,tc.Opts).build('stages',{'prelim'});
            T.feature_b(:) = 10; writetable(T,fullfile(tc.Folder,'metadata_test.csv'));
            [data,~] = INIT([tc.Folder filesep],obj.opts,obj.model);
            tc.verifyEqual(data.X,obj.model.data.Xraw);
            tc.verifyEqual(size(data.X,2),3);
        end
        function testSparseMissingDataRejected(tc)
            file = fullfile(tc.Folder,'metadata.csv'); T = readtable(file);
            T.feature_a(1) = NaN; writetable(T,file);
            tc.verifyError(@() InstanceSpace(tc.Folder,tc.Opts).build('stages',{'prelim'}),'ISA:INIT:incompleteData');
            T.feature_a(1) = 1; T.algo_a(1) = NaN; writetable(T,file);
            tc.verifyError(@() InstanceSpace(tc.Folder,tc.Opts).build('stages',{'prelim'}),'ISA:INIT:incompleteData');
        end
        function testSiftedRestoresInput(tc)
            file = fullfile(tc.Folder,'metadata.csv'); T = readtable(file);
            T.algo_a = T.feature_a; T.algo_b = T.feature_c; writetable(T,file);
            opts = tc.Opts; opts.auto.preproc = false; opts.perf.AbsPerf = true; opts.perf.epsilon = .5;
            opts.sifted.flag = true; opts.sifted.rho = .99;
            obj = InstanceSpace(tc.Folder,opts).build('stages',{'prelim','sifted'});
            tc.verifyEqual(size(obj.model.data.X,2),2);
            obj.opts.sifted.flag = false;
            obj = obj.build('stages',{'sifted'});
            opts.sifted.flag = false;
            fresh = InstanceSpace(tc.Folder,opts).build('stages',{'prelim','sifted'});
            tc.verifyEqual(obj.model.data,fresh.model.data);
            tc.verifyEqual(obj.model.featsel.idx,1:4);
        end
        function testPartialSaveResume(tc)
            obj = InstanceSpace(tc.Folder,tc.Opts);
            obj.opts.pythia.skip = true;
            stages = {'init','prelim','sifted','pilot','cloister','pythia','trace'};
            for i = 1:numel(stages)
                obj = obj.build('stages',stages(i)); obj.save();
                loaded = InstanceSpace.load(tc.Folder);
                tc.verifyEqual(loaded.completedStages,obj.completedStages);
                tc.verifyEqual(loaded.model.data,obj.model.data);
                tc.verifyEqual(loaded.model.opts,obj.model.opts);
                obj = loaded;
            end
        end
        function testImplicitInitCallback(tc)
            stages = {}; snapshots = {};
            obj = InstanceSpace(tc.Folder,tc.Opts).build('stages',{'prelim'},'onStage',@capture);
            tc.verifyEqual(stages,{'init','prelim'});
            tc.verifyEqual(snapshots{1}.data,obj.model.init.data);
            tc.verifyFalse(isfield(snapshots{1},'prelim'));
            tc.verifyTrue(isfield(snapshots{2},'prelim'));
            function capture(stage,snapshot)
                stages{end+1} = stage;
                snapshots{end+1} = snapshot;
            end
        end
        function testInitStageAndPrelimReplay(tc)
            obj = InstanceSpace(tc.Folder,tc.Opts).build('stages',{'init'});
            tc.verifyEqual(obj.completedStages,{'init'});
            tc.verifyEqual(obj.model.data,obj.model.init.data);
            tc.verifyFalse(isfield(obj.model,'prelim'));
            input = obj.model.init.data;
            obj.save(); loaded = InstanceSpace.load(tc.Folder);
            tc.verifyEqual(loaded.completedStages,{'init'});
            tc.verifyEqual(loaded.model.init.data,input);
            obj = obj.build('stages',{'prelim','sifted','pilot'});
            expected = obj.model.data;
            % PRELIM must neither reload the file nor transform its own output.
            delete(fullfile(tc.Folder,'metadata.csv'));
            replay = obj.build('stages',{'prelim'});
            tc.verifyEqual(replay.model.data,expected);
            tc.verifyEqual(replay.model.init.data,input);
            tc.verifyEqual(replay.completedStages,{'init','prelim'});
            tc.verifyFalse(isfield(replay.model,'pilot'));
            resumed = loaded.build('stages',{'prelim'});
            tc.verifyEqual(resumed.model.data,expected);
        end
        function testInitReloadInvalidatesAndTracksOptions(tc)
            obj = InstanceSpace(tc.Folder,tc.Opts).build('stages',{'prelim','sifted','pilot'});
            obj.opts.prelim.nanThreshold = .9;
            tc.verifyError(@() obj.build('stages',{'prelim'}),'ISA:InstanceSpace:staleOptions');
            obj = obj.build('stages',{'init'});
            tc.verifyEqual(obj.completedStages,{'init'});
            tc.verifyFalse(isfield(obj.model,'prelim'));
            tc.verifyFalse(isfield(obj.model,'pilot'));
            tc.verifyFalse(isfield(obj.model,'featsel'));
            file = fullfile(tc.Folder,'metadata.csv'); T = readtable(file);
            T.feature_a = T.feature_a + 5; writetable(T,file);
            previous = obj.model.init.data.X;
            obj = obj.build('stages',{'init'});
            tc.verifyEqual(obj.model.init.data.X(:,1),previous(:,1)+5,'AbsTol',1e-12);
            obj.opts.selvars.feats = {'feature_a','feature_b'};
            tc.verifyError(@() obj.build('stages',{'prelim'}),'ISA:InstanceSpace:staleOptions');
            obj = obj.build('stages',{'init','prelim'});
            tc.verifyEqual(size(obj.model.init.data.X,2),2);
        end
        function testLegacyModelWithoutInitSnapshot(tc)
            opts = tc.Opts; opts.pythia.skip = true;
            obj = InstanceSpace(tc.Folder,opts).build();
            expected = obj.model.data;
            obj.model = rmfield(obj.model,'init');
            obj.model.stageOptions = rmfield(obj.model.stageOptions,'init');
            obj.save();
            copyfile(fullfile(tc.Folder,'metadata.csv'),fullfile(tc.Folder,'metadata_test.csv'));
            loaded = InstanceSpace.load(tc.Folder);
            tc.verifyEqual(loaded.completedStages,obj.completedStages);
            % Old models can still explore without original training metadata.
            movefile(fullfile(tc.Folder,'metadata.csv'),fullfile(tc.Folder,'training.csv'));
            loaded = loaded.explore(tc.Folder);
            tc.verifyEqual(numel(loaded.testResults),1);
            movefile(fullfile(tc.Folder,'training.csv'),fullfile(tc.Folder,'metadata.csv'));
            loaded = loaded.build('stages',{'prelim'});
            tc.verifyEqual(loaded.model.data,expected);
            tc.verifyTrue(isfield(loaded.model,'init'));
        end
        function testStaleOptions(tc)
            obj = InstanceSpace(tc.Folder,tc.Opts).build('stages',{'prelim','sifted','pilot'});
            obj.opts.norm.flag = false;
            tc.verifyError(@() obj.build('stages',{'pilot'}),'ISA:InstanceSpace:staleOptions');
            obj = obj.build('stages',{'prelim','sifted','pilot'});
            tc.verifyFalse(obj.model.opts.norm.flag);
            obj.opts.pilot.dims = 7;
            tc.verifyError(@() obj.build('stages',{'pilot'}),'ISA:ISAvalidateOpts:notMember');
        end
    end
end
