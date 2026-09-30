function issue63_prepare_inputs(root, outputFile)
% Prepare fixed PILOT inputs with the reference exporter's seed and options.
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
work=tempname; mkdir(work); cleanup=onCleanup(@() rmdir(work,'s')); %#ok<NASGU>
copyfile(fullfile(root,'test','data','metadata.csv'),fullfile(work,'metadata.csv'));
opts=struct('general',struct('seed',42,'parallel',false,'verbose',false), ...
    'outputs',struct('csv',false,'png',false,'fig',false,'web',false));
obj=InstanceSpace(work,opts).build('stages',{'prelim','sifted'});
X=obj.model.data.X; Y=obj.model.data.Y; labels=obj.model.data.featlabels; pilotOpts=obj.opts.pilot;
save(outputFile,'X','Y','labels','pilotOpts');
fprintf('ISSUE63_INPUTS: %d instances, %d features, %d algorithms\n',size(X,1),size(X,2),size(Y,2));
end
