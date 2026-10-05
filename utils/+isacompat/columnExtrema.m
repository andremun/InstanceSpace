function [lo, hi] = columnExtrema(X)
% columnExtrema  Column minima/maxima with NaNs omitted.
% Scoped to nonempty real double matrices used by PRELIM and CLOISTER.
% An all-NaN column returns NaN; Inf values are retained, not treated as missing.
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
    lo = min(X, [], 1, 'omitnan');
    if nargout > 1, hi = max(X, [], 1, 'omitnan'); end
    return
end
validateattributes(X, {'double'}, {'real','2d','nonempty'});
lo = NaN(1,size(X,2));
hi = lo;
for j = 1:size(X,2)
    values = X(~isnan(X(:,j)),j);
    if ~isempty(values)
        lo(j) = min(values);
        hi(j) = max(values);
    end
end
end
