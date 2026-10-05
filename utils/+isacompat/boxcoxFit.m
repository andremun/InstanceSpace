function [transformed, lambda] = boxcoxFit(x)
% boxcoxFit  Fit the Gaussian profile-likelihood Box-Cox power transform.
% MATLAB retains Financial Toolbox boxcox. Octave fits the published profile
% likelihood with fminsearch; a constant sample uses lambda=1 (unidentifiable).
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

if ~isacompat.isOctave()
    [transformed, lambda] = boxcox(x);
    return
end
validateattributes(x, {'double'}, {'real','vector','nonempty','finite','positive'});
if all(x == x(1))
    lambda = 1;
else
    logs = log(x(:));
    centered = logs-mean(logs);
    objective = @(power) profileLoss(power, centered);
    [lambda,~,flag] = fminsearch(objective,0, ...
        optimset('Display','off','TolX',1e-8,'TolFun',1e-8,'MaxIter',1000));
    if flag <= 0
        error('ISA:compat:boxcoxConvergence','Box-Cox profile-likelihood optimization did not converge.');
    end
end
transformed = isacompat.boxcoxApply(x,lambda);
end

function loss = profileLoss(lambda, centered)
% Centering log(x) removes the Jacobian term algebraically and reduces
% overflow. It does not change the maximizing lambda (scale invariance).
if abs(lambda)<1e-10
    values = centered;
else
    values = expm1(lambda*centered)/lambda;
end
variance = mean((values-mean(values)).^2);
if ~isfinite(variance) || variance<=0
    loss = Inf;
else
    loss = numel(centered)/2*log(variance);
end
end
