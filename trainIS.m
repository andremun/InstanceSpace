function [model] = trainIS(rootdir)
% -------------------------------------------------------------------------
% trainIS.m
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
scriptdisc('trainIS.m');
% -------------------------------------------------------------------------
% Collect all the data from the files
disp(['Root Directory: ' rootdir]);
datafile = [rootdir 'metadata.csv'];
optsfile = [rootdir 'options.json'];
if ~isfile(datafile) || ~isfile(optsfile)
    error(['Please place the datafiles in the directory ''' rootdir '''']);
end
opts = jsondecode(fileread(optsfile));
disp('-------------------------------------------------------------------------');
disp('-> Listing options to be used:');
optfields = fieldnames(opts);
for i = 1:length(optfields)
    disp(optfields{i});
    disp(opts.(optfields{i}));
end
disp('-------------------------------------------------------------------------');
disp('-> Loading the data.');
Xbar = readtable(datafile);
varlabels = Xbar.Properties.VariableNames;
isname = strcmpi(varlabels,'instances');
isfeat = strncmpi(varlabels,'feature_',8);
isalgo = strncmpi(varlabels,'algo_',5);
issource = strcmpi(varlabels,'source');
instlabels = Xbar{:,isname};
if isnumeric(instlabels)
    instlabels = num2cell(instlabels);
    instlabels = cellfun(@(x) num2str(x),instlabels,'UniformOutput',false); %#ok<NASGU>
end
if any(issource)
    S = categorical(Xbar{:,issource}); %#ok<NASGU>
end
X = Xbar{:,isfeat};
Y = Xbar{:,isalgo};
% -------------------------------------------------------------------------
% PROBABLY HERE SHOULD DO A SANITY CHECK, I.E., IS THERE TOO MANY NANS?
idx = all(isnan(X),2) | all(isnan(Y),2);
X = X(~idx,:);
Y = Y(~idx,:);
idx = mean(isnan(X),1)>=0.20; % These features are very weak.
if any(idx)
    warning('-> There are features with too many missing values. They are being removed to increase speed.');
    X = X(:,~idx);
end
ninst = size(X,1);
nuinst = size(unique(X,'rows'),1);
if nuinst/ninst<0.5
    warning('-> There are too many repeated instances. It is unlikely that this run will produce good results.');
end
% -------------------------------------------------------------------------
% Giving the oportunity to pick and choose which features/algorithms to
% work with
featlabels = varlabels(isfeat);
if isfield(opts,'selvars') && isfield(opts.selvars,'feats')
    disp('-------------------------------------------------------------------------');
    msg = '-> Using the following features: ';
    isselfeat = false(1,length(featlabels));
    for i=1:length(opts.selvars.feats)
        isselfeat = isselfeat | strcmp(featlabels,opts.selvars.feats{i});
        msg = [msg opts.selvars.feats{i} ' ']; %#ok<AGROW>
    end
    disp(msg);
    X = X(:,isselfeat);
    featlabels = featlabels(isselfeat);
end

algolabels = varlabels(isalgo);
if isfield(opts,'selvars') && isfield(opts.selvars,'algos')
    disp('-------------------------------------------------------------------------');
    msg = '-> Using the following algorithms: ';
    isselalgo = false(1,length(algolabels));
    for i=1:length(opts.selvars.algos)
        isselalgo = isselalgo | strcmp(algolabels,opts.selvars.algos{i});
        msg = [msg opts.selvars.algos{i} ' ']; %#ok<AGROW>
    end
    disp(msg);
    Y = Y(:,isselalgo);
    algolabels = algolabels(isselalgo);
end
nalgos = size(Y,2);
% -------------------------------------------------------------------------
% Storing the raw data for further processing, e.g., graphs
Xraw = X; %#ok<NASGU>
Yraw = Y;
% -------------------------------------------------------------------------
% Removing the template data such that it can be used in the labels of
% graphs and figures.
featlabels = strrep(featlabels,'feature_','');
algolabels = strrep(algolabels,'algo_','');
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
if opts.perf.MaxPerf
    [rankPerf,rankAlgo] = sort(Y,2,'descend');
    bestPerformace = rankPerf(:,1);
    P = rankAlgo(:,1);
    if opts.perf.AbsPerf
        Ybin = Y>=opts.perf.epsilon;
        msg = [msg 'higher than ' num2str(opts.perf.epsilon)];
    else
        Ybin = bsxfun(@ge,Y,(1-opts.perf.epsilon).*bestPerformace); % One is good, zero is bad
        msg = [msg 'within ' num2str(round(100.*opts.perf.epsilon)) '% of the best.'];
    end
else
    [rankPerf,rankAlgo] = sort(Y,2,'ascend');
    bestPerformace = rankPerf(:,1);
    P = rankAlgo(:,1);
    if opts.perf.AbsPerf
        Ybin = Y<=opts.perf.epsilon;
        msg = [msg 'less than ' num2str(opts.perf.epsilon)];
    else
        Ybin = bsxfun(@le,Y,(1+opts.perf.epsilon).*bestPerformace);
        msg = [msg 'within ' num2str(round(100.*opts.perf.epsilon)) '% of the best.'];
    end
end
W = abs(Y-bestPerformace);
W(W==0) = min(W(W~=0));
disp(msg);

idx = all(Ybin==0,1);
if any(idx)
    warning('-> There are algorithms with no ''good'' instances. They are being removed to increase speed.');
    Yraw = Yraw(:,~idx);
    Y = Y(:,~idx);
    Ybin = Ybin(:,~idx);
    algolabels = algolabels(~idx);
    nalgos = size(Y,2);
end
% -------------------------------------------------------------------------
% Testing for ties. If there is a tie in performance, we pick an algorithm
% at random.
bestAlgos = bsxfun(@eq,Y,bestPerformace);
multipleBestAlgos = sum(bestAlgos,2)>1;
aidx = 1:nalgos;
for i=1:size(Y,1)
    if multipleBestAlgos(i)
        aux = aidx(bestAlgos(i,:));
        P(i) = aux(randi(length(aux),1)); % Pick one at random
    end
end
disp(['-> For ' num2str(round(100.*mean(multipleBestAlgos))) '% of the instances there is ' ...
      'more than one best algorithm. Random selection is used to break ties.']);
numGoodAlgos = sum(Ybin,2);
beta = numGoodAlgos>opts.general.betaThreshold*nalgos;
% -------------------------------------------------------------------------
% Automated pre-processing
if opts.auto.preproc
    disp('=========================================================================');
    disp('-> Auto-pre-processing.');
    disp('=========================================================================');
    % Eliminate extreme outliers, i.e., any point that exceedes 5 times the
    % inter quantile range, by bounding them to that value.
    if opts.bound.flag
        [X, model.bound] = boundOutliers(X);
    end
    % Normalize the data using Box-Cox and Z-transformations
    if opts.norm.flag
        [X, Y, model.norm] = autoNormalize(X, Y);
    end
end
% -------------------------------------------------------------------------
% If we are only meant to take some observations
disp('-------------------------------------------------------------------------');
ninst = size(X,1);
fractional = isfield(opts,'selvars') && ...
             isfield(opts.selvars,'smallscaleflag') && ...
             opts.selvars.smallscaleflag && ...
             isfield(opts.selvars,'smallscale') && ...
             isfloat(opts.selvars.smallscale);
fileindexed = isfield(opts,'selvars') && ...
              isfield(opts.selvars,'fileidxflag') && ...
              opts.selvars.fileidxflag && ...
              isfield(opts.selvars,'fileidx') && ...
              isfile(opts.selvars.fileidx);
if fractional
    disp(['-> Creating a small scale experiment for validation. Percentage of subset: ' ...
        num2str(round(100.*opts.selvars.smallscale,2)) '%']);
    aux = cvpartition(ninst,'HoldOut',opts.selvars.smallscale);
    subsetIndex = aux.test;
elseif fileindexed
    disp('-> Using a subset of the instances.');
    subsetIndex = false(size(X,1),1);
    aux = table2array(readtable(opts.selvars.fileidx));
    aux(aux>ninst) = [];
    subsetIndex(aux) = true;
else
    disp('-> Using the complete set of the instances.');
    subsetIndex = true(ninst,1);
end

if fileindexed || fractional
    X = X(subsetIndex,:);
    Y = Y(subsetIndex,:);
    Ybin = Ybin(subsetIndex,:);
    beta = beta(subsetIndex);
    bestPerformace = bestPerformace(subsetIndex); 
    P = P(subsetIndex);
    W = W(subsetIndex,:);
end
nfeats = size(X,2);
% -------------------------------------------------------------------------
% Automated feature selection.
% Keep track of the features that have been removed so we can use them
% later
model.featsel.idx = 1:nfeats;
if opts.auto.featsel
    disp('=========================================================================');
    disp('-> Auto-feature selection.');
    disp('=========================================================================');
    % Detect correlations between features and algorithms. Keep the top
    % CORTHRESHOLD correlated features for each algorithm
    if opts.corr.flag
        [X, model.corr] = checkCorrelation(X, Y, opts.corr);
        featlabels = featlabels(model.corr.selvars);
        model.featsel.idx = model.featsel.idx(model.corr.selvars);
    end
    % Detect similar features, by clustering them according to their
    % correlation. We assume that the lowest value possible is best, as
    % this will improve the projection into two dimensions. We set a hard
    % limit of 10 features. The selection criteria is an average silhouete
    % value above 0.65
    if opts.clust.flag
        [X, model.clust] = clusterFeatureSelection(X, Ybin, opts.clust);
        featlabels = featlabels(model.clust.selvars);
        model.featsel.idx = model.featsel.idx(model.clust.selvars);
    end
end
% -------------------------------------------------------------------------
% This is the final subset of features. Calculate the two dimensional
% projection using the PILOT algorithm (Munoz et al. Mach Learn 2018)
disp('=========================================================================');
disp('-> Calling PILOT to find the optimal projection.');
disp('=========================================================================');
model.pilot = PILOT(X, Y, featlabels, opts.pbldr);
% -------------------------------------------------------------------------
% Finding the empirical bounds based on the ranges of the features and the
% correlations of the different edges.
disp('=========================================================================');
disp('-> Finding empirical bounds using CLOISTER.');
disp('=========================================================================');
model.cloist = CLOISTER(X, model.pilot.A, opts.sbound);
% -------------------------------------------------------------------------
% Algorithm selection. Fit a model that would separate the space into
% classes of good and bad performance. 
disp('=========================================================================');
disp('-> Summoning PYTHIA to train the prediction models.');
disp('=========================================================================');
model.pythia = PYTHIA(model.pilot.Z, Yraw(subsetIndex,:), Ybin, W, bestPerformace, algolabels, opts.oracle);
% -------------------------------------------------------------------------
% Calculating the algorithm footprints.
disp('=========================================================================');
disp('-> Calling TRACE to perform the footprint analysis.');
disp('=========================================================================');
if opts.footprint.usesim
    disp('  -> TRACE will use PYTHIA''s results to calculate the footprints.');
    model.trace = TRACE(model.pilot.Z, model.pythia.Yhat, model.pythia.selection0, beta, algolabels, opts.footprint);
else
    disp('  -> TRACE will use experimental data to calculate the footprints.');
    model.trace = TRACE(model.pilot.Z, Ybin, P, beta, algolabels, opts.footprint);
end
% -------------------------------------------------------------------------
% Preparing the outputs for further analysis
scriptfcn; % Loading all the subfunctions into memory

if opts.outputs.csv
    % Storing the output data as a CSV files. This is for easier
    % post-processing. All workspace data will be stored in a matlab file
    % later.
    scriptcsv;
    if opts.outputs.web
        scriptweb;
    end
end
% -------------------------------------------------------------------------
% Making all the plots. First, plotting the features and performance as
% scatter plots.
if opts.outputs.png
    scriptpng;
end
% -------------------------------------------------------------------------
model.opts = opts;
disp('-------------------------------------------------------------------------');
disp('-> Storing the raw MATLAB results for post-processing and/or debugging.');
save([rootdir 'model.mat'],'-struct','model'); % Save the main results
save([rootdir 'workspace.mat']); % Save the full workspace for debugging
disp(['-> Completed! Elapsed time: ' num2str(toc(startProcess)) 's']);
disp('EOF:SUCCESS');
end
