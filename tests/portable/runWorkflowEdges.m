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
function runWorkflowEdges()
% Optional subset routes and archive failure/round-trip contracts.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
opts=struct('general',struct('verbose',false),'outputs',struct('csv',false,'png',false), ...
    'selvars',struct('smallscaleflag',true,'smallscale',.8));
a=InstanceSpace(fullfile(root,'test','data'),opts); a=a.build('stages',{'prelim'});
assert(size(a.model.data.X,1)>100 && size(a.model.data.X,1)<212);
opts.selvars=struct('densityflag',true,'mindistance',.1);
a=InstanceSpace(fullfile(root,'test','data'),opts); a=a.build('stages',{'prelim'});
assert(all(isfinite(a.model.data.X(:))));
if isacompat.isOctave()
    filename=[tempname '.mat']; cleanup=onCleanup(@() delete(filename));
    original=struct('values',[1 NaN 3]); isacompat.saveModel(filename,original);
    expectError(@() isacompat.saveModel(filename,struct('bad',@sin)),'ISA:compat:archiveType');
    assert(isequaln(isacompat.loadModel(filename),original));
    expectError(@() isacompat.requireFeature('parallel'),'ISA:compat:unsupportedFeature');
end
fprintf('[PORTABLE] PASS: fractional/density subsets and archive failure safety.\n');
end
function expectError(f,id)
try, f(); catch err, assert(strcmp(err.identifier,id)); return; end
error('ISA:portable:missingError','Expected %s.',id);
end
