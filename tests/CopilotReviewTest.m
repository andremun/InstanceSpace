classdef CopilotReviewTest < matlab.unittest.TestCase
% CopilotReviewTest  Regression cases for PR 64 review findings.
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
        Dims = {2, 3};
        Boundary = {'populated', 'empty', 'absent'};
    end
    methods (Test)
        function testLegacyPLSBoundaryMigration(tc, Dims, Boundary)
            state = rng; tc.addTeardown(@() rng(state)); rng(27);
            X = rand(30,4) + [2 4 6 8]; Y = rand(30,2);
            opts = ISAdefaults(struct()); opts.pilot.method = 'pls'; opts.pilot.dims = Dims;
            opts.pilot.verbose = false;
            model.opts = opts;
            model.data.X = X;
            model.pilot = PILOT(X,Y,{'a','b','c','d'},opts.pilot);
            fitted = model.pilot;
            model.pilot = rmfield(model.pilot,'Xmean');
            if strcmp(Boundary,'populated')
                model.cloist = CLOISTER(X,model.pilot.A,opts.cloister);
            elseif strcmp(Boundary,'empty')
                model.cloist = struct('Zedge',[],'Zecorr',zeros(0,Dims));
            end
            migrated = ISAmigrateModel(model);
            tc.verifyEqual(migrated.pilot.Xmean,fitted.Xmean);
            tc.verifyEqual(migrated.pilot.Z,fitted.Z);
            tc.verifyEqual(migrated.pilot.A,fitted.A);
            if strcmp(Boundary,'populated')
                shift = fitted.Xmean*fitted.A';
                for field = {'Zedge','Zecorr'}
                    tc.verifyEqual(migrated.cloist.(field{1}),model.cloist.(field{1})-shift,'AbsTol',1e-12);
                end
                tc.verifyEqual(migrated.cloist.ZedgeFaces,model.cloist.ZedgeFaces);
                tc.verifyEqual(migrated.cloist.ZecorrFaces,model.cloist.ZecorrFaces);
            elseif strcmp(Boundary,'empty')
                tc.verifyEqual(migrated.cloist,model.cloist);
            else
                tc.verifyFalse(isfield(migrated,'cloist'));
            end
            tc.verifyEqual(ISAmigrateModel(migrated),migrated, ...
                'Loading an already migrated model must not translate vertices again.');
        end
        function testSerialSiftedDoesNotDispatchToExistingPool(tc)
            pool = gcp('nocreate');
            if isempty(pool)
                pool = parpool('local',2,'SpmdEnabled',false);
                tc.addTeardown(@() delete(pool));
            end
            state = rng; tc.addTeardown(@() rng(state));
            opts = ISAdefaults(struct()); opts.sifted.parallel = false;
            % This path resets the cache but has too few features for GA.
            X = [(1:10)',(10:-1:1)',ones(10,1)]; Y = (1:10)';
            profile on;
            stopProfile = onCleanup(@() profile('off')); %#ok<NASGU>
            [selected,out] = SIFTED(X,Y,Y<5,{'a','b','c'},opts.sifted);
            info = profile('info');
            tc.verifyFalse(any(contains({info.FunctionTable.FunctionName},'parfevalOnAll')), ...
                'Serial SIFTED must not queue worker cache resets.');
            tc.verifyEqual(selected,X); tc.verifyEqual(out.selvars,1:3);
            tc.verifyEqual(gcp('nocreate'),pool);
        end
    end
end
