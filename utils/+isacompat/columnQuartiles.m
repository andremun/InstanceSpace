function [med, spread] = columnQuartiles(X)
% columnQuartiles  NaN-omitting median and IQR along the first dimension.
% MATLAB uses its native operations. Octave explicitly selects quantile
% method 5 (linear interpolation at sample midpoints), independent of its
% default method. All-missing columns yield NaN; constant columns have IQR 0.
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
    med = nanmedian(X, 1);
    spread = iqr(X, 1);
    return
end
validateattributes(X, {'double'}, {'real','2d','nonempty'});
med = NaN(1,size(X,2));
spread = med;
for j = 1:size(X,2)
    values = X(~isnan(X(:,j)),j);
    if ~isempty(values)
        q = quantile(values, [0.25 0.5 0.75], 1, 5);
        med(j) = q(2);
        spread(j) = q(3)-q(1);
    end
end
end
