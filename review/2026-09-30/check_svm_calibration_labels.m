function differences = check_svm_calibration_labels()
% Compare raw and calibrated SVM decisions before considering calibration removal.
% Returns the number of changed training predictions for a fixed synthetic case.
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
%
%   Munoz, M.A., Villanova, L., Baatar, D. & Smith-Miles, K. (2018).
%   Instance spaces for machine learning classification. Machine
%   Learning, 107(1), 109-147. https://doi.org/10.1007/s10994-017-5629-5
% -------------------------------------------------------------------------

state = rng;
cleanup = onCleanup(@() rng(state)); %#ok<NASGU>
rng(7);
Z = rand(40,2); y = rand(40,1) > .5;
raw = fitcsvm(Z,y,'KernelFunction','gaussian');
calibrated = fitSVMPosterior(raw);
differences = nnz(predict(raw,Z) ~= predict(calibrated,Z));
fprintf('CALIBRATED_LABEL_DIFFERENCES=%d\n',differences);
end
