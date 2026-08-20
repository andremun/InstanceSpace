function root = testRepoRoot()
% testRepoRoot  Absolute path (trailing slash) to the repository root, one
% level up from tests/. Used instead of a relative './test/data/' literal
% because MATLAB's current working directory during actual matlab.unittest
% execution is not guaranteed to be the directory test_integration.m was
% launched from -- observed directly in CI (#39): TestClassSetup methods
% ran with cwd set to tests/ itself, turning './test/data/...' into
% 'tests/test/data/...' and failing every copyfile/readtable call.
% mfilename('fullpath') always refers to THIS file's own location
% regardless of who calls it, so this is robust to cwd wherever it's
% called from.

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

root = fileparts(fileparts(mfilename('fullpath')));
if ~(endsWith(root, '/') || endsWith(root, '\'))
    root = [root filesep];
end
end
