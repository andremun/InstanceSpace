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
        function testSignedRelativePerformance(tc)
            p = struct('MaxPerf',false,'AbsPerf',false,'epsilon',0.05,'betaThreshold',0.5,'auto',false);
            X = (1:4)'; Y = [-10 -1; 0 0; 0 1; -1 1];
            [~,loss,out] = PRELIM(X,Y,p);
            tc.verifyEqual(out.Ybin,logical([1 0;1 1;1 0;1 0]));
            tc.verifyEqual(out.Ybest,[-10;0;0;-1]);
            tc.verifyTrue(all(loss(:)>=0));
            p.MaxPerf = true;
            [~,other,maxout] = PRELIM(X,-Y,p);
            tc.verifyEqual(other,loss); tc.verifyEqual(maxout.Ybin,out.Ybin);
            tc.verifyEqual(maxout.Ybest,-out.Ybest);
        end
    end
end
