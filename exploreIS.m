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
fprintf('[EXPLORE] Root directory: %s\n', rootdir);
modelfile = [rootdir 'model.mat'];
datafile = [rootdir 'metadata_test.csv'];
if ~isfile(modelfile) || ~isfile(datafile)
    error(['Please place the datafiles in the directory ''' rootdir '''']);
end
model = load(modelfile);
model = ISAmigrateModel(model);        % migrate opts.oracle->pythia, svm/knn->classifiers, etc.
model.opts = ISAdefaults(model.opts);  % fill in any defaults absent from the saved model
if model.opts.general.verbose
    fprintf('[EXPLORE] Listing options in use:\n');
    optfields = fieldnames(model.opts);
    for i = 1:length(optfields)
        fprintf('%s\n', optfields{i});
        disp(model.opts.(optfields{i}));
    end
end
fprintf('[EXPLORE] Loading metadata_test.csv.\n');
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
    aux = find(strcmpi(strtrim(out.data.algolabels{ii}), strtrim(model.data.algolabels)));
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
fprintf('[EXPLORE] Calculating the binary measure of performance.\n');
msg = 'An algorithm is good if its performance is ';
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
    out.data.Ybest = rankPerf(:,1);
    out.data.P = rankAlgo(:,1);
    if model.opts.perf.AbsPerf
        out.data.Ybin = out.data.Y>=model.opts.perf.epsilon;
        msg = [msg 'higher than ' num2str(model.opts.perf.epsilon)];
    else
        out.data.Ybest(out.data.Ybest==0) = eps;
        out.data.Y(out.data.Y==0) = eps;
        out.data.Y = 1-bsxfun(@rdivide,out.data.Y,out.data.Ybest);
        out.data.Ybin = (1-bsxfun(@rdivide,Yaux,out.data.Ybest))<=model.opts.perf.epsilon;
        msg = [msg 'within ' num2str(round(100.*model.opts.perf.epsilon)) '% of the best.'];
    end
else
    Yaux = out.data.Y;
    Yaux(isnan(Yaux)) = Inf;
    [rankPerf,rankAlgo] = sort(Yaux,2,'ascend');
    out.data.Ybest = rankPerf(:,1);
    out.data.P = rankAlgo(:,1);
    if model.opts.perf.AbsPerf
        out.data.Ybin = out.data.Y<=model.opts.perf.epsilon;
        msg = [msg 'less than ' num2str(model.opts.perf.epsilon)];
    else
        out.data.Ybest(out.data.Ybest==0) = eps;
        out.data.Y(out.data.Y==0) = eps;
        out.data.Y = bsxfun(@rdivide,out.data.Y,out.data.Ybest)-1;
        out.data.Ybin = (bsxfun(@rdivide,Yaux,out.data.Ybest)-1)<=model.opts.perf.epsilon;
        msg = [msg 'within ' num2str(round(100.*model.opts.perf.epsilon)) '% of the best.'];
    end
end
fprintf('[EXPLORE] %s\n', msg);
out.data.numGoodAlgos = sum(out.data.Ybin,2);
out.data.beta = out.data.numGoodAlgos>model.opts.perf.betaThreshold*nalgos;
% ---------------------------------------------------------------------
% Automated pre-processing
if model.opts.auto.preproc && model.opts.bound.flag
    fprintf('[EXPLORE] Auto-pre-processing. Bounding outliers, scaling and normalizing the data.\n');
    fprintf('[EXPLORE] Removing extreme outliers from the feature values.\n');
    himask = bsxfun(@gt, out.data.X, model.prelim.hibound);
    lomask = bsxfun(@lt, out.data.X, model.prelim.lobound);
    out.data.X = out.data.X.*~(himask | lomask) + bsxfun(@times, himask, model.prelim.hibound) + ...
                                                  bsxfun(@times, lomask, model.prelim.lobound);
end

if model.opts.auto.preproc && model.opts.norm.flag
    fprintf('[EXPLORE] Auto-normalizing the data.\n');
    out.data.X = bsxfun(@minus, out.data.X, model.prelim.minX) + 1;
    out.data.X(~isnan(out.data.X) & out.data.X < 1) = 1;  % clamp to training minimum
    for ii = 1:length(model.prelim.lambdaX)
        lambda = model.prelim.lambdaX(ii);
        x = out.data.X(:,ii);
        idx = ~isnan(x);
        if abs(lambda) < 1e-10
            % lambda ≈ 0: log transform
            x(idx) = log(x(idx));
        else
            x(idx) = (x(idx).^lambda - 1) ./ lambda;
        end
        out.data.X(:,ii) = x;
    end
    out.data.X = bsxfun(@rdivide, bsxfun(@minus, out.data.X, model.prelim.muX), model.prelim.sigmaX);

    % If the algorithm is new, something else should be made...
    out.data.Y = (out.data.Y - model.prelim.minY) + eps;
    out.data.Y(out.data.Y <= 0) = eps;  % clamp any remaining non-positives
    for ii = 1:modelalgos
        lambda = model.prelim.lambdaY(ii);
        y = out.data.Y(:,ii);
        idx = ~isnan(y);
        if abs(lambda) < 1e-10
            y(idx) = log(y(idx));
        else
            y(idx) = (y(idx).^lambda - 1) ./ lambda;
        end
        out.data.Y(:,ii) = y;
    end
    out.data.Y(:,1:modelalgos) = bsxfun(@rdivide, bsxfun(@minus, out.data.Y(:,1:modelalgos), model.prelim.muY), model.prelim.sigmaY);
    if newalgos>0
        [~,out.data.Y(:,modelalgos+1:nalgos),out.norm] = autoNormalize(ones(ninst,1), ... % Dummy variable
                                                                       out.data.Y(:,modelalgos+1:nalgos));
    end
end
if ~isreal(out.data.X)
    error('ISA:exploreIS:complexX', ...
        'Feature matrix X is complex after normalisation. Check test data range vs training data.');
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
% model.opts.pythia is guaranteed here: ISAmigrateModel (line 25) has already
% renamed opts.oracle->opts.pythia for pre-v1.7 models, and ISAdefaults (line 26)
% fills in any missing pythia sub-fields, so this call is always safe.
out.pythia = PYTHIA(out.pilot.Z, out.data.Yraw, out.data.Ybin, out.data.Ybest, ...
                    out.data.algolabels, model.opts.pythia, model.pythia);
% -------------------------------------------------------------------------
% Validating the footprints (evaluation mode: polygons from training reused)
out.trace = TRACE(out.pilot.Z, out.data.Ybin, out.pythia.Yhat, out.data.P, ...
                  out.data.beta, out.data.algolabels, model.opts.trace, model.trace);

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

fprintf('[EXPLORE] Storing the raw MATLAB results for post-processing and/or debugging.\n');
save([rootdir 'workspace_test.mat']); % Save the full workspace for debugging
fprintf('[EXPLORE] Completed in %.1f s.\n', toc(startProcess));
fprintf('EOF:SUCCESS\n');
end

% =========================================================================
%  SUBFUNCTIONS
% =========================================================================

function [X, Y, out] = autoNormalize(X, Y)
% autoNormalize  Fit and apply Box-Cox + Z-score normalisation to X and Y.
%
% Used only for algorithms present in the test metadata but absent from
% the training model (out.data.Y columns modelalgos+1:nalgos): those
% columns have no pre-fit model.prelim.lambdaY/muY/sigmaY to reuse, so a
% fresh normalisation is fit directly on the test data for them.
nfeats = size(X, 2);
nalgos = size(Y, 2);
out.minX = min(X, [], 1, 'omitnan');
X = bsxfun(@minus, X, out.minX) + 1;
out.lambdaX = zeros(1, nfeats);
out.muX = zeros(1, nfeats);
out.sigmaX = zeros(1, nfeats);
for i = 1:nfeats
    aux = X(:,i);
    idx = isnan(aux);
    [aux, out.lambdaX(i)] = boxcox(aux(~idx));
    [aux, out.muX(i), out.sigmaX(i)] = zscore(aux);
    X(~idx,i) = aux;
end

out.minY = min(Y(:), [], 'omitnan');
Y = (Y - out.minY) + eps;
out.lambdaY = zeros(1, nalgos);
out.muY = zeros(1, nalgos);
out.sigmaY = zeros(1, nalgos);
for i = 1:nalgos
    aux = Y(:,i);
    idx = isnan(aux);
    [aux, out.lambdaY(i)] = boxcox(aux(~idx));
    [aux, out.muY(i), out.sigmaY(i)] = zscore(aux);
    Y(~idx,i) = aux;
end
end
% =========================================================================
