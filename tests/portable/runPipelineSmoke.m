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
function runPipelineSmoke(dims,withPlots)
% End-to-end release contract; isolated deterministic metadata and outputs.
if nargin<1, dims=2; end
if nargin<2, withPlots=false; end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
folder=fullfile(root,'test','data','octave-validation',sprintf('pipeline%d',dims));
if ~isfolder(folder), mkdir(folder); end
for testing=[false true]
    if testing, name='metadata_test.csv'; else, name='metadata.csv'; end
    fid=fopen(fullfile(folder,name),'w');
    fprintf(fid,'instances,feature_a,feature_b,feature_c,feature_d,feature_e,feature_f,algo_left,algo_right,source');
    if testing, fprintf(fid,',algo_extra'); end
    fprintf(fid,'\n');
    for k=1:48
        t=k/8+double(testing)/32;
        x=[t,sin(t),cos(t),sin(2*t),cos(3*t),t.^2];
        y=[2+sin(t),2-sin(t)];
        fprintf(fid,'"instância,%d",%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,group%d',k,x,y,mod(k,2));
        if testing, fprintf(fid,',%.15g',3+cos(t)); end
        fprintf(fid,'\n');
    end
    fclose(fid);
end
opts=struct('general',struct('verbose',false), ...
    'sifted',struct('K',4,'Replicates',2,'diagnostics',false), ...
    'pilot',struct('dims',dims,'analytic',true), ...
    'pythia',struct('classifier','knn','tuning','sobol','nTuningIter',4,'kFold',3), ...
    'outputs',struct('csv',true,'png',withPlots,'fig',false));
obj=InstanceSpace(folder,opts); obj=obj.build();
assert(size(obj.model.pilot.Z,2)==dims);
assert(strcmp(obj.model.data.instlabels{1},'instância,1'));
assert(isfinite(obj.model.trace.space.measure) && obj.model.trace.space.measure>0);
loaded=InstanceSpace.load(folder);
assert(isequal(loaded.model.pythia.Yhat,obj.model.pythia.Yhat));
obj=obj.explore(folder); loaded=loaded.explore(folder);
a=obj.testResults{end}; b=loaded.testResults{end};
assert(size(a.data.Y,2)==3 && all(isfinite(a.data.Y(:))));
assert(isequal(a.pythia.Yhat,b.pythia.Yhat));
assert(max(abs(a.pythia.Pr0hat(:)-b.pythia.Pr0hat(:)))<1e-12);
assert(isequaln(a.trace.summary,b.trace.summary));
assert(isequal(cellstr(loaded.model.data.S),cellstr(obj.model.data.S)));
assert(numel(dir(fullfile(folder,'*.csv')))>5);
if withPlots
    plots=dir(fullfile(folder,'*.png')); assert(numel(plots)>5);
    fid=fopen(fullfile(folder,plots(1).name),'r'); signature=fread(fid,8,'uint8')'; fclose(fid);
    assert(isequal(signature,[137 80 78 71 13 10 26 10]));
end
fprintf('[PORTABLE] PASS: %dD full build/save/load/explore, CSV and PNG=%d.\n',dims,withPlots);
end
