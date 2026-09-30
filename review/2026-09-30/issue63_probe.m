function issue63_probe(root, inputFile, outputPrefix)
% Run two fixed-input PILOT solves and save candidates and execution metadata.
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



addpath(root,fullfile(root,'core'),fullfile(root,'utils'),fullfile(root,'output'));
in=load(inputFile); opts=in.pilotOpts; opts.verbose=false;
meta=struct('release',version,'platform',computer,'threads',maxNumCompThreads, ...
    'blas',version('-blas'),'ntries',opts.ntries,'seed',opts.seed,'root',root);
assert(isempty(gcp('nocreate')),'Probe requires no parallel pool.');
for repeat=1:2
    out=PILOT(in.X,in.Y,in.labels,opts);
    [~,winner]=max(out.perf); sorted=sort(out.perf,'descend'); gap=sorted(1)-sorted(2);
    name=sprintf('%s-%d.mat',outputPrefix,repeat);
    save(name,'out','meta','winner','gap');
    fprintf('ISSUE63 %s repeat %d: threads=%d winner=%d gap=%.17g\n',outputPrefix,repeat,meta.threads,winner,gap);
end
end
