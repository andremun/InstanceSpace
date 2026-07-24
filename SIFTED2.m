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
%   Simpson, C., Munoz, M.A., Kandanaarachchi, S. & Campello, R.J.G.B.
%   (2025). ISA3: A 3-dimensional expansion of Instance Space Analysis.
%   Machine Learning, 114, 240. https://doi.org/10.1007/s10994-025-06871-5
%
%   Smith-Miles, K. & Munoz, M.A. (2023). Instance Space Analysis for
%   Algorithm Testing. ACM Computing Surveys, 55(12), Article 255.
%   https://doi.org/10.1145/3572895
% -------------------------------------------------------------------------
function [X, out] = SIFTED2(X, Y, Ybin, featlabels, opts)
% SIFTED2  Deprecated. Use SIFTED instead.
%
% This wrapper exists for backward compatibility. All new code should call
% SIFTED directly; SIFTED2 was promoted and renamed as part of the v1.7
% refactor (spec Phase 6).
warning('ISA:SIFTED2:deprecated', ...
    'SIFTED2 is deprecated and will be removed in a future release. Call SIFTED instead.');
[X, out] = SIFTED(X, Y, Ybin, featlabels, opts);
end
