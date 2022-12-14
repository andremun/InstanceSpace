function model = buildIS(rootdir)
% -------------------------------------------------------------------------
% buildIS.m
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
scriptdisc('buildIS.m');
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
useparallel = isfield(opts,'parallel') && isfield(opts.parallel,'flag') && opts.parallel.flag;
if useparallel
    disp('-------------------------------------------------------------------------');
    disp('-> Starting parallel processing pool.');
    delete(gcp('nocreate'));
    if  isfield(opts.parallel,'ncores') && isnumeric(opts.parallel.ncores)
        mypool = parpool('local',opts.parallel.ncores,'SpmdEnabled',false);
    else
        mypool = parpool('local','SpmdEnabled',false);
    end
    if ispc
        addAttachedFiles(mypool,{'svmpredict.mexw64','svmtrain.mexw64'});
    elseif isunix
        addAttachedFiles(mypool,{'svmpredict.mexa64','svmtrain.mexa64'});
    elseif ismac
        addAttachedFiles(mypool,{'libsvmpredict.mexmaci64','libsvmtrain.mexmaci64'});
    end
end
disp('-------------------------------------------------------------------------');
disp('-> Loading the data.');
Xbar = readtable(datafile);
varlabels = Xbar.Properties.VariableNames;
isname = strcmpi(varlabels,'instances');
isfeat = strncmpi(varlabels,'feature_',8);
isalgo = strncmpi(varlabels,'algo_',5);
issource = strcmpi(varlabels,'source');
model.data.instlabels = Xbar{:,isname};
if isnumeric(model.data.instlabels)
    model.data.instlabels = num2cell(model.data.instlabels);
    model.data.instlabels = cellfun(@(x) num2str(x),model.data.instlabels,'UniformOutput',false);
end
if any(issource)
    model.data.S = categorical(Xbar{:,issource});
end
model.data.X = Xbar{:,isfeat};
model.data.Y = Xbar{:,isalgo};
% -------------------------------------------------------------------------
% Giving the oportunity to pick and choose which features/algorithms to
% work with
model.data.featlabels = varlabels(isfeat);
if isfield(opts,'selvars') && isfield(opts.selvars,'feats')
    disp('-------------------------------------------------------------------------');
    msg = '-> Using the following features: ';
    isselfeat = false(1,length(model.data.featlabels));
    for i=1:length(opts.selvars.feats)
        isselfeat = isselfeat | strcmp(model.data.featlabels,opts.selvars.feats{i});
        msg = [msg opts.selvars.feats{i} ' ']; %#ok<AGROW>
    end
    disp(msg);
    model.data.X = model.data.X(:,isselfeat);
    model.data.featlabels = model.data.featlabels(isselfeat);
end

model.data.algolabels = varlabels(isalgo);
if isfield(opts,'selvars') && isfield(opts.selvars,'algos')
    disp('-------------------------------------------------------------------------');
    msg = '-> Using the following algorithms: ';
    isselalgo = false(1,length(model.data.algolabels));
    for i=1:length(opts.selvars.algos)
        isselalgo = isselalgo | strcmp(model.data.algolabels,opts.selvars.algos{i});
        msg = [msg opts.selvars.algos{i} ' ']; %#ok<AGROW>
    end
    disp(msg);
    model.data.Y = model.data.Y(:,isselalgo);
    model.data.algolabels = model.data.algolabels(isselalgo);
end
% nalgos = size(model.data.Y,2);
% -------------------------------------------------------------------------
% PROBABLY HERE SHOULD DO A SANITY CHECK, I.E., IS THERE TOO MANY NANS?
idx = all(isnan(model.data.X),2) | all(isnan(model.data.Y),2);
if any(idx)
    warning('-> There are instances with too many missing values. They are being removed to increase speed.');
    model.data.X = model.data.X(~idx,:);
    model.data.Y = model.data.Y(~idx,:);
    model.data.instlabels = model.data.instlabels(~idx);
    if isfield(model.data,'S')
        model.data.S = model.data.S(~idx);
    end
end
idx = mean(isnan(model.data.X),1)>=0.20; % These features are very weak.
if any(idx)
    warning('-> There are features with too many missing values. They are being removed to increase speed.');
    model.data.X = model.data.X(:,~idx);
    model.data.featlabels = model.data.featlabels(~idx);
end
ninst = size(model.data.X,1);
nuinst = size(unique(model.data.X,'rows'),1);
if nuinst/ninst<0.5
    warning('-> There are too many repeated instances. It is unlikely that this run will produce good results.');
end
% -------------------------------------------------------------------------
% Storing the raw data for further processing, e.g., graphs
model.data.Xraw = model.data.X;
model.data.Yraw = model.data.Y;
% -------------------------------------------------------------------------
% Removing the template data such that it can be used in the labels of
% graphs and figures.
model.data.featlabels = strrep(model.data.featlabels,'feature_','');
model.data.algolabels = strrep(model.data.algolabels,'algo_','');
% -------------------------------------------------------------------------
% Running PRELIM as to pre-process the data, including scaling and bounding
opts.prelim = opts.perf;
opts.prelim.bound = opts.bound.flag;
opts.prelim.norm = opts.norm.flag;
[model.data.X,model.data.Y,model.data.Ybest,...
    model.data.Ybin,model.data.P,model.data.numGoodAlgos,...
    model.data.beta,model.prelim] = PRELIM(model.data.X, model.data.Y, opts.prelim);

idx = all(~model.data.Ybin,1);
if any(idx)
    warning('-> There are algorithms with no ''good'' instances. They are being removed to increase speed.');
    model.data.Yraw = model.data.Yraw(:,~idx);
    model.data.Y = model.data.Y(:,~idx);
    model.data.Ybin = model.data.Ybin(:,~idx);
    model.data.algolabels = model.data.algolabels(~idx);
    nalgos = size(model.data.Y,2);
    if nalgos==0
        error('-> There are no ''good'' algorithms. Please verify the binary performance measure. STOPPING!')
    end
end
% -------------------------------------------------------------------------
% If we are only meant to take some observations
disp('-------------------------------------------------------------------------');
ninst = size(model.data.X,1);
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
bydensity = isfield(opts,'selvars') && ...
            isfield(opts.selvars,'densityflag') && ...
            opts.selvars.densityflag && ...
            isfield(opts.selvars,'mindistance') && ...
            isfloat(opts.selvars.mindistance) && ...
            isfield(opts.selvars,'type') && ...
            ischar(opts.selvars.type);
if fractional
    disp(['-> Creating a small scale experiment for validation. Percentage of subset: ' ...
        num2str(round(100.*opts.selvars.smallscale,2)) '%']);
    state = rng;
    rng('default');
    aux = cvpartition(ninst,'HoldOut',opts.selvars.smallscale);
    rng(state);
    subsetIndex = aux.test;
elseif fileindexed
    disp('-> Using a subset of the instances.');
    subsetIndex = false(size(model.data.X,1),1);
    aux = table2array(readtable(opts.selvars.fileidx));
    aux(aux>ninst) = [];
    subsetIndex(aux) = true;
elseif bydensity
    disp('-> Creating a small scale experiment for validation based on density.');
    subsetIndex = FILTER(model.data.X, model.data.Y, model.data.Ybin, ...
                         opts.selvars);
    subsetIndex = ~subsetIndex;
    disp(['-> Percentage of instances retained: ' ...
          num2str(round(100.*mean(subsetIndex),2)) '%']);
else
    disp('-> Using the complete set of the instances.');
    subsetIndex = true(ninst,1);
end

if fileindexed || fractional || bydensity
    if bydensity
        model.data_dense = model.data;
    end
    model.data.X = model.data.X(subsetIndex,:);
    model.data.Y = model.data.Y(subsetIndex,:);
    model.data.Xraw = model.data.Xraw(subsetIndex,:);
    model.data.Yraw = model.data.Yraw(subsetIndex,:);
    model.data.Ybin = model.data.Ybin(subsetIndex,:);
    model.data.beta = model.data.beta(subsetIndex);
    model.data.numGoodAlgos = model.data.numGoodAlgos(subsetIndex);
    model.data.Ybest = model.data.Ybest(subsetIndex); 
    model.data.P = model.data.P(subsetIndex);
    model.data.instlabels = model.data.instlabels(subsetIndex);
    if isfield(model.data,'S')
        model.data.S = model.data.S(subsetIndex);
    end
end
nfeats = size(model.data.X,2);
% -------------------------------------------------------------------------
% Automated feature selection.
% Keep track of the features that have been removed so we can use them
% later
model.featsel.idx = 1:nfeats;
if opts.sifted.flag
    disp('=========================================================================');
    disp('-> Calling SIFTED for auto-feature selection.');
    disp('=========================================================================');
    [model.data.X, model.sifted] = SIFTED(model.data.X, model.data.Y, model.data.Ybin, opts.sifted);
    model.data.featlabels = model.data.featlabels(model.sifted.selvars);
    model.featsel.idx = model.featsel.idx(model.sifted.selvars);

    if bydensity
        disp('-> Creating a small scale experiment for validation based on density.');
        % model.data.featlabels = model.data_dense.featlabels(model.sifted.selvars);
        subsetIndex = FILTER(model.data_dense.X(:,model.featsel.idx), ...
                             model.data_dense.Y, ...
                             model.data_dense.Ybin, ...
                             opts.selvars);
        subsetIndex = ~subsetIndex;
        model.data.X = model.data_dense.X(subsetIndex,model.featsel.idx);
        model.data.Y = model.data_dense.Y(subsetIndex,:);
        model.data.Xraw = model.data_dense.Xraw(subsetIndex,:);
        model.data.Yraw = model.data_dense.Yraw(subsetIndex,:);
        model.data.Ybin = model.data_dense.Ybin(subsetIndex,:);
        model.data.beta = model.data_dense.beta(subsetIndex);
        model.data.numGoodAlgos = model.data_dense.numGoodAlgos(subsetIndex);
        model.data.Ybest = model.data_dense.Ybest(subsetIndex);
        model.data.P = model.data_dense.P(subsetIndex);
        model.data.instlabels = model.data_dense.instlabels(subsetIndex);
        if isfield(model.data_dense,'S')
            model.data.S = model.data_dense.S(subsetIndex);
        end
        disp(['-> Percentage of instances retained: ' ...
              num2str(round(100.*mean(subsetIndex),2)) '%']);
    end
end
% -------------------------------------------------------------------------
% This is the final subset of features. Calculate the two dimensional
% projection using the PILOT algorithm (Munoz et al. Mach Learn 2018)
disp('=========================================================================');
disp('-> Calling PILOT to find the optimal projection.');
disp('=========================================================================');
model.pilot = PILOT(model.data.X, model.data.Y, model.data.featlabels, opts.pilot);
% -------------------------------------------------------------------------
% Finding the empirical bounds based on the ranges of the features and the
% correlations of the different edges.
disp('=========================================================================');
disp('-> Finding empirical bounds using CLOISTER.');
disp('=========================================================================');
model.cloist = CLOISTER(model.data.X, model.pilot.A, opts.cloister);
% -------------------------------------------------------------------------
% Algorithm selection. Fit a model that would separate the space into
% classes of good and bad performance. 
disp('=========================================================================');
disp('-> Summoning PYTHIA to train the prediction models.');
disp('=========================================================================');
model.pythia = PYTHIA(model.pilot.Z, model.data.Yraw, model.data.Ybin, model.data.Ybest, model.data.algolabels, opts.pythia);
% -------------------------------------------------------------------------
% Calculating the algorithm footprints.
disp('=========================================================================');
disp('-> Calling TRACE to perform the footprint analysis.');
disp('=========================================================================');
if opts.trace.usesim
    disp('  -> TRACE will use PYTHIA''s results to calculate the footprints.');
    model.trace = TRACE(model.pilot.Z, model.pythia.Yhat, model.pythia.selection0, model.data.beta, model.data.algolabels, opts.trace);
else
    disp('  -> TRACE will use experimental data to calculate the footprints.');
    model.trace = TRACE(model.pilot.Z, model.data.Ybin, model.data.P, model.data.beta, model.data.algolabels, opts.trace);
end

if useparallel
    disp('-------------------------------------------------------------------------');
    disp('-> Closing parallel processing pool.');
    delete(mypool);
end
% -------------------------------------------------------------------------
% Preparing the outputs for further analysis
model.opts = opts;
% -------------------------------------------------------------------------
disp('-------------------------------------------------------------------------');
disp('-> Storing the raw MATLAB results for post-processing and/or debugging.');
save([rootdir 'model.mat'],'-struct','model'); % Save the main results
save([rootdir 'workspace.mat']); % Save the full workspace for debugging
% -------------------------------------------------------------------------
if opts.outputs.csv
    % Storing the output data as a CSV files. This is for easier
    % post-processing. All workspace data will be stored in a matlab file
    % later.
    scriptcsv(model,rootdir);
    if opts.outputs.web
        scriptweb(model,rootdir);
    end
end
% -------------------------------------------------------------------------
% Making all the plots. First, plotting the features and performance as
% scatter plots.
if opts.outputs.png
    scriptpng(model,rootdir);
end
% -------------------------------------------------------------------------
disp(['-> Completed! Elapsed time: ' num2str(toc(startProcess)) 's']);
disp('EOF:SUCCESS');
end
