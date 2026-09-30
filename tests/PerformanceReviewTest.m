classdef PerformanceReviewTest < matlab.unittest.TestCase
% PerformanceReviewTest  Regression cases for the September 2026 review.
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
        function testBatchedCornerHull(tc)
            rng(31); X = rand(40,13)-.5; A = rand(3,13);
            opts = ISAdefaults(struct()); opts.cloister.corrThreshold = 1;
            out = CLOISTER(X,A,opts.cloister);
            lo = min(X); hi = max(X);
            ids = (0:2^13-1)'; corners = zeros(numel(ids),13);
            for j = 1:13
                bits = bitget(uint32(ids),j);
                corners(:,j) = lo(j); corners(bits==1,j) = hi(j);
            end
            Z = corners*A'; [~,expected] = convhull(Z);
            [~,actual] = convhull(out.Zedge);
            tc.verifyEqual(actual,expected,'RelTol',1e-10);
            for direction = [eye(3),rand(3,10)]
                tc.verifyEqual(max(out.Zedge*direction),max(Z*direction),'AbsTol',1e-10);
            end
        end
        function testFilterGreedyOrder(tc)
            opts = struct('mindistance',.2,'type','Ftr');
            X = [0; .1; .25; .3; 1; 1]; Y = ones(6,1);
            removed = FILTER(X,Y,true(size(Y)),opts);
            tc.verifyEqual(removed,logical([0;1;0;1;0;1]));
        end
    end
end
