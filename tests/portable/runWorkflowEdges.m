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
    checkStagedOutputRejection(root);
end
fprintf('[PORTABLE] PASS: fractional/density subsets and archive failure safety.\n');
end
function expectError(f,id)
try, f(); catch err, assert(strcmp(err.identifier,id)); return; end
error('ISA:portable:missingError','Expected %s.',id);
end

function checkStagedOutputRejection(root)
folder=tempname; mkdir(folder); cleanup=onCleanup(@() rmdir(folder,'s'));
copyfile(fullfile(root,'test','data','metadata.csv'),fullfile(folder,'metadata.csv'));
opts=struct('general',struct('verbose',false),'outputs',struct('csv',false,'png',false), ...
    'sifted',struct('flag',false),'pilot',struct('analytic',true,'dims',3,'ntries',1), ...
    'pythia',struct('skip',true),'trace',struct('minInstances',1000));
obj=InstanceSpace(folder,opts);
obj=obj.build('stages',{'prelim','sifted','pilot','cloister','pythia'});
assert(~ismember('trace',obj.completedStages));
names={'model.mat','bounds.csv','footprint_previous.png','footprint_previous.fig'};
for j=1:numel(names)
    fid=fopen(fullfile(folder,names{j}),'w'); fprintf(fid,'preserve this existing output'); fclose(fid);
end
for mode={'web','fig'}
    bad=obj; bad.opts.outputs.csv=true;
    if strcmp(mode{1},'web')
        bad.opts.outputs.web=true;
    else
        bad.opts.outputs.png=true; bad.opts.outputs.fig=true;
    end
    before=rng;
    expectError(@() bad.build('stages',{'trace'}),'ISA:compat:unsupportedFeature');
    assert(isequal(rng,before));
    % Even a non-final stage must reject these unsupported output options.
    expectError(@() bad.build('stages',{'prelim'}),'ISA:compat:unsupportedFeature');
    if strcmp(mode{1},'fig')
        container=bad.model; container.opts=bad.opts;
        expectError(@() scriptpng(container,[folder filesep]),'ISA:compat:unsupportedFeature');
    end
    for j=1:numel(names)
        assert(strcmp(fileread(fullfile(folder,names{j})),'preserve this existing output'));
    end
end
% Supported output settings still allow the final stage to complete and save.
obj=obj.build('stages',{'trace'});
assert(ismember('trace',obj.completedStages));
assert(isfield(isacompat.loadModel(fullfile(folder,'model.mat')),'trace'));
checkExploreOutputRejection(obj,folder,names);
fprintf('[PORTABLE] PASS: staged output rejection preserves existing files.\n');
end

function checkExploreOutputRejection(obj,folder,names)
copyfile(fullfile(folder,'metadata.csv'),fullfile(folder,'metadata_test.csv'));
for mode={'web','fig'}
    saved=obj.model; saved.opts.outputs.csv=true;
    if strcmp(mode{1},'web')
        saved.opts.outputs.web=true;
    else
        saved.opts.outputs.png=true; saved.opts.outputs.fig=true;
    end
    % Exercise options retained by an archive, not merely current obj.opts.
    isacompat.saveModel(fullfile(folder,'model.mat'),saved);
    loaded=InstanceSpace.load(folder);
    loaded.opts.outputs.web=false; loaded.opts.outputs.fig=false;
    for j=2:numel(names)
        fid=fopen(fullfile(folder,names{j}),'w'); fprintf(fid,'preserve this existing output'); fclose(fid);
    end
    before=rng;
    expectError(@() loaded.explore(folder,'onStage',@unexpectedStage), ...
        'ISA:compat:unsupportedFeature');
    assert(isequal(rng,before));
    assert(isempty(loaded.testDirs) && isempty(loaded.testResults));
    for j=2:numel(names)
        assert(strcmp(fileread(fullfile(folder,names{j})),'preserve this existing output'));
    end
end
% Training-only options must not stop evaluation of an already trained model.
obj.model.opts.pythia.skip=false;
obj.model.opts.pythia.tuning='bayes';
obj.opts.outputs.web=true; % Explore uses frozen model options, not this edit.
obj=obj.explore(folder);
assert(numel(obj.testResults)==1 && numel(obj.testDirs)==1);
fprintf('[PORTABLE] PASS: explore checks frozen outputs before evaluation or mutation.\n');
end
function unexpectedStage(varargin)
error('ISA:portable:unexpectedStage','Explore evaluated a stage before rejecting output options.');
end
