function data = ISAsubsetData(data, subsetIndex, featIdx)
% ISAsubsetData  Subset rows (and optionally feature columns) of a data struct.
%   data = ISAsubsetData(data, subsetIndex) subsets all row-indexed fields.
%   data = ISAsubsetData(data, subsetIndex, featIdx) also selects feature
%   columns featIdx from data.X (used in the post-SIFTED density path).

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
%
%   Simpson, C., Munoz, M.A., Kandanaarachchi, S. & Campello, R.J.G.B.
%   (2025). ISA3: A 3-dimensional expansion of Instance Space Analysis.
%   Machine Learning, 114, 240. https://doi.org/10.1007/s10994-025-06871-5
% -------------------------------------------------------------------------
if nargin < 3
    data.X = data.X(subsetIndex, :);
else
    data.X = data.X(subsetIndex, featIdx);
    data.featlabels = data.featlabels(featIdx);
end
data.Y            = data.Y(subsetIndex, :);
data.Xraw         = data.Xraw(subsetIndex, :);
data.Yraw         = data.Yraw(subsetIndex, :);
data.Ybin         = data.Ybin(subsetIndex, :);
data.beta         = data.beta(subsetIndex);
data.numGoodAlgos = data.numGoodAlgos(subsetIndex);
data.Ybest        = data.Ybest(subsetIndex);
data.P            = data.P(subsetIndex);
data.instlabels   = data.instlabels(subsetIndex);
if isfield(data, 'S')
    data.S = data.S(subsetIndex);
end
end
