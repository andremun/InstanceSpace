function out = exploreIS(rootdir)
% exploreIS  Thin backward-compatibility wrapper around InstanceSpace.
%
% Requires model.mat to already exist in rootdir (written by buildIS).
% Preserves the original exploreIS(rootdir) calling convention for
% callers -- notably the MATILDA web platform -- that invoke this entry
% point directly. New code should use InstanceSpace directly:
%
%   obj = InstanceSpace.load(rootdir);
%   obj = obj.explore(rootdir);
%   out  = obj.getResults(1);

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

if ~(endsWith(rootdir, '/') || endsWith(rootdir, '\'))
    rootdir = [rootdir '/'];
end
if ~isfile([rootdir 'model.mat']) || ~isfile([rootdir 'metadata_test.csv'])
    error(['Please place the datafiles in the directory ''' rootdir '''']);
end

obj = InstanceSpace.load(rootdir);
obj = obj.explore(rootdir);
out = obj.testResults{end};

fprintf('EOF:SUCCESS\n');
end
