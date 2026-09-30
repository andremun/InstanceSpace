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
        function testSiftedRestoresInput(tc)
            obj = InstanceSpace(tc.Folder,tc.Opts).build('stages',{'prelim','sifted'});
            original = obj.model.data;
            obj.model.data.X = original.X(:,[1 3]);
            obj.model.data.featlabels = original.featlabels([1 3]);
            obj.model.featsel.idx = [1 3];
            obj = obj.build('stages',{'sifted'});
            tc.verifyEqual(obj.model.data,original);
            tc.verifyEqual(obj.model.featsel.idx,1:4);
        end
        function testPartialSaveResume(tc)
            obj = InstanceSpace(tc.Folder,tc.Opts);
            stages = {'prelim','sifted','pilot','cloister'};
            for i = 1:numel(stages)
                obj = obj.build('stages',stages(i)); obj.save();
                loaded = InstanceSpace.load(tc.Folder);
                tc.verifyEqual(loaded.completedStages,obj.completedStages);
                tc.verifyEqual(loaded.model.data,obj.model.data);
                tc.verifyEqual(loaded.model.opts,obj.model.opts);
                obj = loaded;
            end
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
