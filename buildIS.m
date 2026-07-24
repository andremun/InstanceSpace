function model = buildIS(rootdir)
% buildIS  Thin backward-compatibility wrapper around InstanceSpace (spec §7.5).
%
% Preserves the pre-Phase-7 buildIS(rootdir) calling convention (a plain
% function taking a directory and returning the in-memory model struct)
% for callers -- notably the MATILDA web platform -- that invoke this
% entry point directly. New code should use InstanceSpace directly:
%
%   obj = InstanceSpace(rootdir);
%   obj = obj.build();

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

if ~(endsWith(rootdir, '/') || endsWith(rootdir, '\'))
    rootdir = [rootdir '/'];
end
if ~isfile([rootdir 'metadata.csv']) || ~isfile([rootdir 'options.json'])
    error(['Please place the datafiles in the directory ''' rootdir '''']);
end

obj = InstanceSpace(rootdir);
obj = obj.build();
model = obj.model;

fprintf('EOF:SUCCESS\n');
end
