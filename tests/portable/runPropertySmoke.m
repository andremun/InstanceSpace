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
function runPropertySmoke()
% Property validation must reject invalid assignments without changing state.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
oldPath=path; pathCleanup=onCleanup(@() path(oldPath)); %#ok<NASGU>
addpath(root);
folder=tempname; mkdir(folder); cleanup=onCleanup(@() rmdir(folder,'s')); %#ok<NASGU>
obj=InstanceSpace(folder,[],false);
assert(isstruct(obj.opts) && isscalar(obj.opts));
properties={'rootdir','opts','model','testDirs','testResults','completedStages'};
badValues={{42,['ab';'cd']}, ...
    {42,struct([]),repmat(struct('x',1),1,2)}, ...
    {42,struct([]),repmat(struct('x',1),1,2)}, ...
    {42,cell(2,1),cell(1,1,2)}, ...
    {42,cell(2,1),cell(1,1,2)}, ...
    {42,cell(2,1),cell(1,1,2)}};
for j=1:numel(properties)
    name=properties{j}; before=obj.(name);
    for k=1:numel(badValues{j})
        rejected=false;
        try
            obj.(name)=badValues{j}{k};
        catch err
            assert(strcmp(err.identifier,'ISA:InstanceSpace:invalidProperty'));
            rejected=true;
        end
        assert(rejected && isequaln(obj.(name),before));
    end
end
obj.rootdir=[folder filesep]; obj.opts=struct('custom',17); obj.model=struct('data',[1 NaN]);
assert(obj.opts.custom==17 && isequaln(obj.model.data,[1 NaN]));
for j=4:numel(properties)
    name=properties{j}; obj.(name)={'first','second'};
    assert(isequal(obj.(name),{'first','second'}));
    obj.(name)=cell(1,0); assert(isequal(size(obj.(name)),[1 0]));
end
fprintf('[PORTABLE] PASS: property validation preserves state and accepts supported shapes.\n');
end
