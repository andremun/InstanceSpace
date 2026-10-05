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
function requireOutputSupport(opts)
% Validate frozen output options before either build or explore does work.
if ~isacompat.isOctave(), return; end
if opts.outputs.png && opts.pilot.dims==3 && opts.outputs.fig
    error('ISA:compat:unsupportedFeature','Set outputs.fig=false for Octave PNG output.');
end
if opts.outputs.web
    error('ISA:compat:unsupportedFeature','Web palette export is MATLAB-only; Octave supports CSV and PNG.');
end
end
