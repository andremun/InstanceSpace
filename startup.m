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
% startup.m - Add the ISA Toolkit's function directories to the MATLAB path.
%
% Run this once per MATLAB session (e.g. `run('startup.m')`, or just
% `startup` if the repo root is your current folder) before calling any
% toolkit function directly -- PILOT, SIFTED, PRELIM, CLOISTER, PYTHIA,
% TRACE, FILTER, the output writers, etc. -- without going through
% buildIS, exploreIS, or InstanceSpace. Those three already add the
% subdirectories to the path automatically the first time one is used
% in a session, so example.m and test_integration.m need no separate
% startup.m step.
%
% Safe to run more than once.

root = fileparts(mfilename('fullpath'));
subdirs = {'core', 'output', 'utils', 'deprecated'};
for i = 1:numel(subdirs)
    d = fullfile(root, subdirs{i});
    if isfolder(d)
        addpath(d);
    end
end
clear root subdirs d i

fprintf('[STARTUP] ISA Toolkit function directories added to the MATLAB path.\n');
