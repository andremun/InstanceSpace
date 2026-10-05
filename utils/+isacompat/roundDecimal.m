function y = roundDecimal(x, digits)
% roundDecimal  Round a display matrix to a nonnegative decimal precision.
% MATLAB retains its native implementation. Octave's one-argument round
% rounds ties away from zero. This adapter is for summaries, not computation.
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
    y = round(x, digits);
    return
end
validateattributes(digits, {'numeric'}, {'scalar','integer','>=',0,'<=',15});
scale = 10.^digits;
y = x;
% At magnitudes where scaling would overflow, decimal rounding cannot
% change the represented value. Preserve Inf and NaN as well.
mask = isfinite(x) & abs(x) <= realmax ./ scale;
y(mask) = round(x(mask) .* scale) ./ scale;
end
