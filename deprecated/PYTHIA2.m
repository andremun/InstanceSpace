function out = PYTHIA2(Z, Y, Ybin, Ybest, algolabels, opts)
% PYTHIA2  Deprecated. Use PYTHIA with opts.classifier = 'knn'.
%
% This wrapper exists for backward compatibility. All new code should call
% PYTHIA directly with the desired classifier type set in opts.classifier.

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
warning('ISA:PYTHIA2:deprecated', ...
    ['PYTHIA2 is deprecated and will be removed in a future release. ' ...
     'Call PYTHIA with opts.classifier = ''knn'' instead.']);
opts.classifier = 'knn';
out = PYTHIA(Z, Y, Ybin, Ybest, algolabels, opts);
end
