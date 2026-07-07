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
fprintf('[BUILD] Root directory: %s\n', rootdir);
datafile = [rootdir 'metadata.csv'];
optsfile = [rootdir 'options.json'];
if ~isfile(datafile) || ~isfile(optsfile)
    error(['Please place the datafiles in the directory ''' rootdir '''']);
end
opts = jsondecode(fileread(optsfile));
opts = ISAdefaults(opts);
rng(opts.general.seed, 'twister');
if opts.general.verbose
    fprintf('[BUILD] Listing options in use:\n');
    optfields = fieldnames(opts);
    for i = 1:length(optfields)
        fprintf('%s\n', optfields{i});
        disp(opts.(optfields{i}));
    end
end
if opts.parallel.flag
    fprintf('[BUILD] Starting parallel processing pool.\n');
    delete(gcp('nocreate'));
    if isnumeric(opts.parallel.ncores)
        mypool = parpool('local', opts.parallel.ncores, 'SpmdEnabled', false);
    else
        mypool = parpool('local', 'SpmdEnabled', false);
    end
end
fprintf('[BUILD] Loading metadata.csv.\n');
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
    msg = 'Using the following features: ';
    isselfeat = false(1,length(model.data.featlabels));
    for i=1:length(opts.selvars.feats)
        isselfeat = isselfeat | strcmp(model.data.featlabels,opts.selvars.feats{i});
        msg = [msg opts.selvars.feats{i} ' ']; %#ok<AGROW>
    end
    fprintf('[BUILD] %s\n', msg);
    model.data.X = model.data.X(:,isselfeat);
    model.data.featlabels = model.data.featlabels(isselfeat);
end

model.data.algolabels = varlabels(isalgo);
if isfield(opts,'selvars') && isfield(opts.selvars,'algos')
    msg = 'Using the following algorithms: ';
    isselalgo = false(1,length(model.data.algolabels));
    for i=1:length(opts.selvars.algos)
        isselalgo = isselalgo | strcmp(model.data.algolabels,opts.selvars.algos{i});
        msg = [msg opts.selvars.algos{i} ' ']; %#ok<AGROW>
    end
    fprintf('[BUILD] %s\n', msg);
    model.data.Y = model.data.Y(:,isselalgo);
    model.data.algolabels = model.data.algolabels(isselalgo);
end
% -------------------------------------------------------------------------
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
fprintf('[PRELIM] Calling PRELIM for data pre-processing.\n');
opts.prelim = opts.perf;
opts.prelim.auto = opts.auto.preproc;
opts.prelim.bound = opts.bound.flag;
opts.prelim.norm = opts.norm.flag;
[model.data.X, model.data.Y, model.prelim] = PRELIM(model.data.X, model.data.Y, opts.prelim);
model.data.Ybest        = model.prelim.Ybest;
model.data.Ybin         = model.prelim.Ybin;
model.data.P            = model.prelim.P;
model.data.numGoodAlgos = model.prelim.numGoodAlgos;
model.data.beta         = model.prelim.beta;

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
ninst = size(model.data.X,1);
fractional = opts.selvars.smallscaleflag && isfloat(opts.selvars.smallscale);
fileindexed = opts.selvars.fileidxflag && isfile(opts.selvars.fileidx);
bydensity   = opts.selvars.densityflag && ...
              isfloat(opts.selvars.mindistance) && ...
              ischar(opts.selvars.type);
if fractional
    fprintf('[BUILD] Creating a small scale experiment for validation. Percentage of subset: %s%%\n', ...
        num2str(round(100.*opts.selvars.smallscale,2)));
    state = rng;
    rng('default');
    aux = cvpartition(ninst,'HoldOut',opts.selvars.smallscale);
    rng(state);
    subsetIndex = aux.test;
elseif fileindexed
    fprintf('[BUILD] Using a subset of the instances.\n');
    subsetIndex = false(size(model.data.X,1),1);
    aux = table2array(readtable(opts.selvars.fileidx));
    aux(aux>ninst) = [];
    subsetIndex(aux) = true;
elseif bydensity
    fprintf('[BUILD] Creating a small scale experiment for validation based on density.\n');
    subsetIndex = FILTER(model.data.X, model.data.Y, model.data.Ybin, ...
                         opts.selvars);
    subsetIndex = ~subsetIndex;
    fprintf('[BUILD] Percentage of instances retained: %s%%\n', ...
        num2str(round(100.*mean(subsetIndex),2)));
else
    fprintf('[BUILD] Using the complete set of the instances.\n');
    subsetIndex = true(ninst,1);
end

if fileindexed || fractional || bydensity
    if bydensity
        model.data_dense = model.data;
    end
    model.data = ISAsubsetData(model.data, subsetIndex);
end
nfeats = size(model.data.X,2);
% -------------------------------------------------------------------------
% Automated feature selection.
% Keep track of the features that have been removed so we can use them
% later
model.featsel.idx = 1:nfeats;
if opts.sifted.flag
    fprintf('[SIFTED] Calling SIFTED for automated feature selection.\n');
    [model.data.X, model.sifted] = SIFTED2(model.data.X, model.data.Y, model.data.Ybin, model.data.featlabels, opts.sifted);
    model.data.featlabels = model.data.featlabels(model.sifted.selvars);
    model.featsel.idx = model.featsel.idx(model.sifted.selvars);

    if bydensity
        fprintf('[SIFTED] Creating a small scale experiment for validation based on density.\n');
        subsetIndex = FILTER(model.data_dense.X(:,model.featsel.idx), ...
                             model.data_dense.Y, ...
                             model.data_dense.Ybin, ...
                             opts.selvars);
        subsetIndex = ~subsetIndex;
        model.data = ISAsubsetData(model.data_dense, subsetIndex, model.featsel.idx);
        fprintf('[SIFTED] Percentage of instances retained: %s%%\n', ...
            num2str(round(100.*mean(subsetIndex),2)));
    end
end
% -------------------------------------------------------------------------
% This is the final subset of features. Calculate the two dimensional
% projection using the PILOT algorithm (Munoz et al. Mach Learn 2018)
fprintf('[PILOT] Calling PILOT to find the optimal projection.\n');
model.pilot = PILOT(model.data.X, model.data.Y, model.data.featlabels, opts.pilot);
% -------------------------------------------------------------------------
% Finding the empirical bounds based on the ranges of the features and the
% correlations of the different edges.
fprintf('[CLOISTER] Finding empirical bounds using CLOISTER.\n');
model.cloist = CLOISTER(model.data.X, model.pilot.A, opts.cloister);
% -------------------------------------------------------------------------
% Algorithm selection. Fit a model that would separate the space into
% classes of good and bad performance.
fprintf('[PYTHIA] Summoning PYTHIA to train the prediction models.\n');
model.pythia = PYTHIA(model.pilot.Z, model.data.Yraw, model.data.Ybin, model.data.Ybest, model.data.algolabels, opts.pythia);
% -------------------------------------------------------------------------
% Calculating the algorithm footprints.
fprintf('[TRACE] Calling TRACE to perform the footprint analysis.\n');
opts.trace.pythiaSkip = opts.pythia.skip;  % Yhat is NaN-placeholder shaped, not empty, when skipped
model.trace = TRACE(model.pilot.Z, model.data.Ybin, model.pythia.Yhat, model.data.P, model.data.beta, model.data.algolabels, opts.trace);

if opts.parallel.flag
    fprintf('[BUILD] Closing parallel processing pool.\n');
    delete(mypool);
end
% -------------------------------------------------------------------------
% Preparing the outputs for further analysis
model.opts = opts;
% -------------------------------------------------------------------------
fprintf('[BUILD] Storing the raw MATLAB results for post-processing and/or debugging.\n');
save([rootdir 'model.mat'],'-struct','model'); % Save the main results
save([rootdir 'workspace.mat']); % Save the full workspace for debugging
% -------------------------------------------------------------------------
if opts.outputs.csv
    scriptcsv(model,rootdir);
    if opts.outputs.web
        scriptweb(model,rootdir);
    end
end
% -------------------------------------------------------------------------
if opts.outputs.png
    scriptpng(model,rootdir);
end
% -------------------------------------------------------------------------
fprintf('[BUILD] Completed in %.1f s.\n', toc(startProcess));
fprintf('EOF:SUCCESS\n');
end
