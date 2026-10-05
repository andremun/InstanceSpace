function [x, value] = minimize(objective, initial, viewpoint)
% minimize  Serial quasi-Newton optimizer with explicit runtime option mapping.
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

if nargin<3, viewpoint=false; end
if isacompat.isOctave()
    options = optimset('Display','off');
    if viewpoint
        options = optimset(options,'MaxIter',30000,'TolFun',1e-20);
    end
else
    options = optimoptions('fminunc','Algorithm','quasi-newton','Display','off','UseParallel',false);
    if viewpoint
        options = optimoptions(options,'MaxIterations',30000,'FunctionTolerance',1e-20);
    end
end
[x,value] = fminunc(objective,initial,options);
end
