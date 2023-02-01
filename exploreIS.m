function out = exploreIS(rootdir)
% -------------------------------------------------------------------------
% exploreIS.m
% -------------------------------------------------------------------------
%
% By: Mario Andres Munoz Acosta
%     School of Mathematics and Statistics
%     The University of Melbourne
%     Australia
%     2019
%
% -------------------------------------------------------------------------

startProcess = tic;
scriptdisc('exploreIS.m');
% -------------------------------------------------------------------------
% Collect all the data from the files
disp(['Root Directory: ' rootdir]);
modelfile = [rootdir 'model.mat'];
datafile = [rootdir 'metadata_test.csv'];
if ~isfile(modelfile) || ~isfile(datafile)
    error(['Please place the datafiles in the directory ''' rootdir '''']);
end
model = load(modelfile);
disp('-------------------------------------------------------------------------');
disp('Listing options used:');
optfields = fieldnames(model.opts);
for i = 1:length(optfields)
    disp(optfields{i});
    disp(model.opts.(optfields{i}));
end
disp('-------------------------------------------------------------------------');
disp('-> Loading the data');
Xbar = readtable(datafile);
varlabels = Xbar.Properties.VariableNames;
isname = strcmpi(varlabels,'instances');
isfeat = strncmpi(varlabels,'feature_',8);
isalgo = strncmpi(varlabels,'algo_',5);
issource = strcmpi(varlabels,'source');
out.data.instlabels = Xbar{:,isname};
if isnumeric(out.data.instlabels)
    out.data.instlabels = num2cell(out.data.instlabels);
    out.data.instlabels = cellfun(@(x) num2str(x),out.data.instlabels,'UniformOutput',false);
end
if any(issource)
    out.data.S = categorical(Xbar{:,issource});
end
out.data.X = Xbar{:,isfeat};
out.data.Y = Xbar{:,isalgo};
[ninst,nalgos] = size(out.data.Y);
% -------------------------------------------------------------------------
% HERE CHECK IF THE NUMBER OF ALGORITHMS IS THE SAME AS IN THE MODEL. IF
% NOT, CHECK IF THE NAMES OF THE ALGORITHMS ARE THE SAME, IF NOT, MOVE THE
% DATA IN SUCH WAY THAT THE NON-EXISTING ALGORITHMS ARE MADE NAN AND THE
% NEW ALGORITHMS ARE LAST.
out.data.algolabels = strrep(varlabels(isalgo),'algo_','');
algoexist = zeros(1,nalgos);
for ii=1:nalgos
    aux = find(strcmp(out.data.algolabels{ii},model.data.algolabels));
    if ~isempty(aux)
        algoexist(ii) = aux;
    end
end
newalgos = sum(algoexist==0);
modelalgos = length(model.data.algolabels);
Yaux = NaN+ones(ninst, modelalgos+newalgos);
lblaux = model.data.algolabels;
acc = modelalgos+1;
for ii=1:nalgos
    if algoexist(ii)==0
       Yaux(:,acc) = out.data.Y(:,ii);
       lblaux(:,acc) = out.data.algolabels(ii);
       acc = acc+1;
    else
        Yaux(:,algoexist(ii)) = out.data.Y(:,ii);
        % lblaux(:,acc) = out.data.algolabels(ii);
    end
end
out.data.Y = Yaux;
out.data.algolabels = lblaux;
nalgos = size(out.data.Y,2);
% -------------------------------------------------------------------------
% Storing the raw data for further processing, e.g., graphs
out.data.Xraw = out.data.X;
out.data.Yraw = out.data.Y;
% -------------------------------------------------------------------------
% Determine whether the performance of an algorithm is a cost measure to
% be minimized or a profit measure to be maximized. Moreover, determine
% whether we are using an absolute threshold as good peformance (the
% algorithm has a performance better than the threshold) or a relative
% performance (the algorithm has a performance that is similar that the
% best algorithm minus a percentage).
disp('-------------------------------------------------------------------------');
disp('-> Calculating the binary measure of performance');
msg = '-> An algorithm is good if its performace is ';
MaxPerf = false;
if isfield(model.opts.perf, 'MaxPerf')
    MaxPerf = model.opts.perf.MaxPerf;
elseif  isfield(model.opts.perf, 'MaxMin')
    MaxPerf = model.opts.perf.MaxMin;
else
    warning('Can not find parameter "MaxPerf" in the trained model. We are assuming that performance metric is needed to be minimized.');
end
if MaxPerf
    Yaux = out.data.Y;
    Yaux(isnan(Yaux)) = -Inf;
    [rankPerf,rankAlgo] = sort(Yaux,2,'descend');
    out.data.bestPerformace = rankPerf(:,1);
    out.data.P = rankAlgo(:,1);
    if model.opts.perf.AbsPerf
        out.data.Ybin = out.data.Y>=model.opts.perf.epsilon;
        msg = [msg 'higher than ' num2str(model.opts.perf.epsilon)];
    else
        out.data.bestPerformace(out.data.bestPerformace==0) = eps;
        out.data.Y(out.data.Y==0) = eps;
        out.data.Y = 1-bsxfun(@rdivide,out.data.Y,out.data.bestPerformace);
        out.data.Ybin = (1-bsxfun(@rdivide,Yaux,out.data.bestPerformace))<=model.opts.perf.epsilon;
        msg = [msg 'within ' num2str(round(100.*model.opts.perf.epsilon)) '% of the best.'];
    end
else
    Yaux = out.data.Y;
    Yaux(isnan(Yaux)) = Inf;
    [rankPerf,rankAlgo] = sort(Yaux,2,'ascend');
    out.data.bestPerformace = rankPerf(:,1);
    out.data.P = rankAlgo(:,1);
    if model.opts.perf.AbsPerf
        out.data.Ybin = out.data.Y<=model.opts.perf.epsilon;
        msg = [msg 'less than ' num2str(model.opts.perf.epsilon)];
    else
        out.data.bestPerformace(out.data.bestPerformace==0) = eps;
        out.data.Y(out.data.Y==0) = eps;
        out.data.Y = bsxfun(@rdivide,out.data.Y,out.data.bestPerformace)-1;
        out.data.Ybin = (bsxfun(@rdivide,Yaux,out.data.bestPerformace)-1)<=model.opts.perf.epsilon;
        msg = [msg 'within ' num2str(round(100.*model.opts.perf.epsilon)) '% of the best.'];
    end
end
disp(msg);
out.data.numGoodAlgos = sum(out.data.Ybin,2);
out.data.beta = out.data.numGoodAlgos>model.opts.perf.betaThreshold*nalgos;
% ---------------------------------------------------------------------
% Automated pre-processing
if model.opts.auto.preproc && model.opts.bound.flag
    disp('-------------------------------------------------------------------------');
    disp('-> Auto-pre-processing. Bounding outliers, scaling and normalizing the data.');
    % Eliminate extreme outliers, i.e., any point that exceedes 5 times the
    % inter quantile range, by bounding them to that value.
    disp('-> Removing extreme outliers from the feature values.');
    himask = bsxfun(@gt,out.data.X,model.bound.hibound);
    lomask = bsxfun(@lt,out.data.X,model.bound.lobound);
    out.data.X = out.data.X.*~(himask | lomask) + bsxfun(@times,himask,model.bound.hibound) + ...
                                                  bsxfun(@times,lomask,model.bound.lobound);
end

if model.opts.auto.preproc && model.opts.norm.flag
    % Normalize the data using Box-Cox and out.pilot.Z-transformations
    disp('-> Auto-normalizing the data.');
    out.data.X = bsxfun(@minus,out.data.X,model.norm.minX)+1;
    out.data.X = bsxfun(@rdivide,bsxfun(@power,out.data.X,model.norm.lambdaX)-1,model.norm.lambdaX);
    out.data.X = bsxfun(@rdivide,bsxfun(@minus,out.data.X,model.norm.muX),model.norm.sigmaX);
    
    % If the algorithm is new, something else should be made...
    out.data.Y(out.data.Y==0) = eps; % Assumes that out.data.Y is always positive and higher than 1e-16
    out.data.Y(:,1:modelalgos) = bsxfun(@rdivide,bsxfun(@power,out.data.Y(:,1:modelalgos),model.norm.lambdaY)-1,model.norm.lambdaY);
    out.data.Y(:,1:modelalgos) = bsxfun(@rdivide,bsxfun(@minus,out.data.Y(:,1:modelalgos),model.norm.muY),model.norm.sigmaY);
    if newalgos>0
        [~,out.data.Y(:,modelalgos+1:nalgos),out.norm] = autoNormalize(ones(ninst,1), ... % Dummy variable
                                                                       out.data.Y(:,modelalgos+1:nalgos));
    end
end
% ---------------------------------------------------------------------
% This is the final subset of features.
out.featsel.idx = model.featsel.idx;
out.data.X = out.data.X(:,out.featsel.idx);
out.data.featlabels = strrep(varlabels(isfeat),'feature_','');
out.data.featlabels = out.data.featlabels(model.featsel.idx);
% ---------------------------------------------------------------------
%  Calculate the two dimensional projection using the PBLDR algorithm
%  (Munoz et al. Mach Learn 2018)
out.pilot.Z = out.data.X*model.pilot.A';
% -------------------------------------------------------------------------
% Algorithm selection. Fit a model that would separate the space into
% classes of good and bad performance. 
out.pythia = PYTHIAtest(model.pythia, out.pilot.Z, out.data.Yraw, ...
                        out.data.Ybin, out.data.bestPerformace, ...
                        out.data.algolabels);
% -------------------------------------------------------------------------
% Validating the footprints
if model.opts.trace.usesim
    out.trace = TRACEtest(model.trace, out.pilot.Z, out.pythia.Yhat, ...
                          out.pythia.selection0, out.data.beta, ...
                          out.data.algolabels);
%     out.trace = TRACE(out.pilot.Z, out.pythia.Yhat, out.pythia.selection0, ...
%                       out.data.beta, out.data.algolabels, model.opts.trace);
else
    out.trace = TRACEtest(model.trace, out.pilot.Z, out.data.Ybin, ...
                          out.data.P, out.data.beta, ...
                          out.data.algolabels);
%     out.trace = TRACE(out.pilot.Z, out.data.Ybin, out.data.P, out.data.beta,...
%                       out.data.algolabels, model.opts.trace);
end

out.opts = model.opts;
% -------------------------------------------------------------------------
% Writing the results
if model.opts.outputs.csv
    scriptcsv(out,rootdir);
    if model.opts.outputs.web
        scriptweb(out,rootdir);
    end
end

if model.opts.outputs.png
    scriptpng(out,rootdir);
end

disp('-------------------------------------------------------------------------');
disp('-> Storing the raw MATLAB results for post-processing and/or debugging.');
save([rootdir 'workspace_test.mat']); % Save the full workspace for debugging
disp(['-> Completed! Elapsed time: ' num2str(toc(startProcess)) 's']);
disp('EOF:SUCCESS');
end
% =========================================================================