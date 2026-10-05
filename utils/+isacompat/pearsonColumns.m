function [rho, pval] = pearsonColumns(X)
% pearsonColumns  Pearson column correlations and two-sided significance.
% Implements only corr(X)'s default Pearson/Rows='all'/Tail='both' contract,
% not rank correlations or pairwise deletion. A NaN or constant column yields
% NaN coefficients and p-values. The Octave adapter requires >=3 observations.
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
    [rho, pval] = corr(X);
    return
end
validateattributes(X, {'double'}, {'real','2d','nonempty'});
n = size(X,1);
if n < 3
    error('ISA:compat:insufficientCorrelationRows', ...
        'Octave Pearson significance requires at least three observations.');
end
rho = corr(X);
valid = ~isnan(rho);
% Roundoff can put a valid Pearson coefficient just outside [-1,1].
rho(valid) = max(-1, min(1, rho(valid)));
pval = NaN(size(rho));
% For t = r*sqrt((n-2)/(1-r^2)), the two-sided Student-t tail is
% I_(1-r^2)((n-2)/2,1/2). This avoids a Statistics-package tcdf dependency.
pval(valid) = betainc(max(0, 1-rho(valid).^2), (n-2)/2, 0.5);
% Self-correlation is exactly one for nonconstant, complete columns.
indices = 1:size(rho,1)+1:numel(rho);
indices = indices(valid(indices));
rho(indices) = 1;
pval(indices) = 0;
end
