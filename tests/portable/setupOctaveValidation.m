function report = setupOctaveValidation(installPackages)
% setupOctaveValidation  Explicitly load the isolated pinned Octave packages.
% setupOctaveValidation(true) downloads, verifies and installs missing packages
% under ignored test/data/octave-validation. The default only loads them.
% This developer helper is never called by startup or library algorithms.
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

if nargin < 1, installPackages = false; end
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'utils'));
if ~isacompat.isOctave() || compare_versions(version,'11.1.0','<')
    error('ISA:compat:runtimeVersion','This validation environment requires Octave >=11.1.0.');
end
envdir = fullfile(root,'test','data','octave-validation');
if ~isfolder(envdir)
    if ~installPackages
        error('ISA:compat:missingPackages','Run setupOctaveValidation(true) to provision the pinned packages.');
    end
    mkdir(envdir);
end
pkg('prefix', fullfile(envdir,'packages'));
pkg('local_list', fullfile(envdir,'octave_packages'));
names = {'datatypes','statistics'};
versions = {'1.5.0','2.0.0'};
urls = { ...
    'https://github.com/pr0m1th3as/datatypes/releases/download/release-1.5.0/datatypes-1.5.0.tar.gz', ...
    'https://github.com/gnu-octave/statistics/releases/download/release-2.0.0/statistics-2.0.0.tar.gz'};
checksums = { ...
    '5d4eb5efe22a1388d7f0a4337539fd4d8a79f33d40fb37a23c4da0bd7d3de32b', ...
    'e82c1d6957885ee444adc91e35defd00665ec3c3e1eeae0bd39f89b7b28a6490'};
for k = 1:numel(names)
    installed = pkg('list');
    matches = cellfun(@(p) strcmp(p.name,names{k}) && strcmp(p.version,versions{k}),installed);
    if ~any(matches)
        if ~installPackages
            error('ISA:compat:missingPackages','Missing %s %s; run setupOctaveValidation(true).',names{k},versions{k});
        end
        archive = fullfile(envdir,[names{k} '-' versions{k} '.tar.gz']);
        if ~isfile(archive), urlwrite(urls{k},archive); end
        fid = fopen(archive,'rb');
        if fid < 0, error('ISA:compat:packageRead','Cannot read %s.',archive); end
        guard = onCleanup(@() fclose(fid));
        bytes = fread(fid,Inf,'*uint8');
        clear guard
        if ~strcmp(hash('sha256',char(bytes')),checksums{k})
            error('ISA:compat:packageChecksum','Checksum mismatch for %s; no installation attempted.',archive);
        end
        pkg('install','-local',archive);
    end
    pkg('load',names{k});
end
report = isacompat.diagnostics();
end
