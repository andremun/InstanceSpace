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
function model=loadModel(filename)
model=load(filename);
if isfield(model,'archiveVersion')
    if ~isacompat.isOctave() || ~isnumeric(model.archiveVersion) || ...
            ~isscalar(model.archiveVersion) || ~ismember(model.archiveVersion,[1 2])
        error('ISA:compat:archiveVersion','This portable archive requires the tested Octave environment (schemas 1 and 2).');
    end
    model=isacompat.archiveValue(model.payload,true,model.archiveVersion);
end
end
