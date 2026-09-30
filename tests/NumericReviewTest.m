classdef NumericReviewTest < matlab.unittest.TestCase
% NumericReviewTest  Regression cases for the September 2026 review.
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


    methods (Test)
        function testPositiveRelativeThresholdCompatibility(tc)
            p = struct('MaxPerf',false,'AbsPerf',false,'epsilon',0.05,'betaThreshold',0.5,'auto',false);
            % Include the 5% boundary and positive best values below eps.
            Y = [100 105 106; 1e-20 1.05e-20 2e-20; 3 4 5];
            for maximize = [false true]
                p.MaxPerf = maximize;
                if maximize
                    best = max(Y,[],2); expected = 1 - Y ./ best;
                else
                    best = min(Y,[],2); expected = Y ./ best - 1;
                end
                [~,relative,out] = PRELIM((1:3)',Y,p);
                tc.verifyEqual(relative,expected);
                tc.verifyEqual(out.Ybin,expected <= p.epsilon);
                tc.verifyEqual(out.Ybest,best);
            end
        end
        function testBoundarySeeds(tc)
            opts = ISAdefaults(struct()); opts.pythia.seed = 2^32-1;
            opts.pythia.params = [1 1]; opts.pythia.kFold = 2;
            Z = [(1:8)',[1;3;2;4;6;5;8;7]]; Y = (1:8)';
            a = PYTHIA(Z,Y,Y<5,Y,{'a'},opts.pythia);
            b = PYTHIA(Z,Y,Y<5,Y,{'a'},opts.pythia);
            tc.verifyEqual(a.Ysub,b.Ysub);
            opts.general.seed = 2^32;
            tc.verifyError(@() ISAvalidateOpts(opts),'ISA:ISAvalidateOpts:seedRange');
        end
        function testPerformanceDomainAndZeroTies(tc)
            p = struct('MaxPerf',false,'AbsPerf',false,'epsilon',0.05,'betaThreshold',0.5,'auto',false);
            X = (1:3)'; Y = [0 0;0 1;1 2];
            for maximize = [false true]
                p.MaxPerf = maximize;
                [~,relative,out] = PRELIM(X,Y,p);
                tc.verifyEqual(out.Ybin(1,:),[true true]);
                tc.verifyEqual(relative(1,:),[0 0]);
                tc.verifyEqual(out.Ybest(1),0);
                tc.verifyEqual(out.Ybin,relative <= p.epsilon);
                for absolute = [false true]
                    p.AbsPerf = absolute;
                    tc.verifyError(@() PRELIM(1,[-1 1],p),'ISA:PRELIM:invalidPerformance');
                    tc.verifyError(@() PRELIM(1,[Inf 1],p),'ISA:PRELIM:invalidPerformance');
                    tc.verifyError(@() PRELIM(1,[-1 1],p,out),'ISA:PRELIM:invalidPerformance');
                end
                p.AbsPerf = false;
            end
        end
        function testZeroReferenceWarning(tc)
            p = struct('MaxPerf',true,'AbsPerf',false,'epsilon',0.05,'betaThreshold',0.5,'auto',false);
            tc.verifyWarning(@() PRELIM(1,[0 0],p),'ISA:PRELIM:manyZeroBest');
        end
    end
end
