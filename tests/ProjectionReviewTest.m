classdef ProjectionReviewTest < matlab.unittest.TestCase
% ProjectionReviewTest  Regression cases for the September 2026 review.
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
        function testRankDeficientFallback(tc)
            X = [(1:8)',2*(1:8)']; Y = sin((1:8)');
            out = PILOT(X,Y,{'a','b'},struct('analytic',true,'verbose',false));
            tc.verifyTrue(all(isfinite(out.Z(:))));
            tc.verifySize(out.Z,[8 2]);
        end
        function testAnalyticAvoidsDistances(tc)
            rng(3); X = rand(30,4); Y = rand(30,2);
            profile on;
            cleanup = onCleanup(@() profile('off'));
            out = PILOT(X,Y,{'a','b','c','d'},struct('analytic',true,'verbose',false));
            info = profile('info');
            names = {info.FunctionTable.FunctionName};
            tc.verifyFalse(any(strcmp(names,'pdist')));
            tc.verifyEqual(out.Z,X*out.A','AbsTol',1e-10);
        end
    end
end
