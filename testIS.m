function [out] = testIS(rootdir)
% -------------------------------------------------------------------------
% testIS.m
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
printdisclaim('testIS.m');
% -------------------------------------------------------------------------
% Collect all the data from the files
modelfile = [rootdir 'results.mat'];
datafile = [rootdir 'metadata.csv'];
optsfile = [rootdir 'options.json'];
if ~isfile(modelfile) || ~isfile(datafile) || ~isfile(optsfile)
    error('Please place the datafiles in the directory specified.');
end
opts = jsondecode(fileread(optsfile));
disp('-------------------------------------------------------------------------');
disp('Listing options used:');
optfields = fieldnames(opts);
for i = 1:length(optfields)
    disp(optfields{i});
    disp(opts.(optfields{i}));
end

disp('-------------------------------------------------------------------------');
disp('-> Loading the data');
model = load(modelfile);
Xbar = readtable(datafile);
varlabels = Xbar.Properties.VariableNames;
isname = strcmpi(varlabels,'instances');
isfeat = strncmpi(varlabels,'feature_',8);
isalgo = strncmpi(varlabels,'algo_',5);
issource = strcmpi(varlabels,'source');
instlabels = Xbar{:,isname};
if isnumeric(instlabels)
    instlabels = num2cell(instlabels);
    instlabels = cellfun(@(x) num2str(x),instlabels,'UniformOutput',false);
end
if any(issource)
    S = categorical(Xbar{:,issource});
    sourcelabels = cellstr(unique(S));
end
X = Xbar{:,isfeat};
Y = Xbar{:,isalgo};
nalgos = size(Y,2);
ninst = size(X,1);
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
if model.opts.perf.MaxMin
    Y(isnan(Y)) = -Inf;
    [bestPerformace,portfolio] = max(Y,[],2);
    if model.opts.perf.AbsPerf
        Ybin = Y>=model.opts.perf.epsilon;
        msg = [msg 'higher than ' num2str(model.opts.perf.epsilon)];
    else
        Ybin = bsxfun(@ge,Y,(1-model.opts.perf.epsilon).*bestPerformace); % One is good, zero is bad
        msg = [msg 'within ' num2str(round(100.*model.opts.perf.epsilon)) '% of the best.'];
    end
else
    Y(isnan(Y)) = Inf;
    [bestPerformace,portfolio] = min(Y,[],2);
    if model.opts.perf.AbsPerf
        Ybin = Y<=model.opts.perf.epsilon;
        msg = [msg 'less than ' num2str(model.opts.perf.epsilon)];
    else
        Ybin = bsxfun(@le,Y,(1+model.opts.perf.epsilon).*bestPerformace);
        msg = [msg 'within ' num2str(round(100.*model.opts.perf.epsilon)) '% of the best.'];
    end
end
disp(msg);
beta = sum(Ybin,2)>model.opts.general.betaThreshold*nalgos;

% ---------------------------------------------------------------------
% Automated pre-processing
if model.opts.bound.flag
    disp('-------------------------------------------------------------------------');
    disp('-> Auto-pre-processing. Bounding outliers, scaling and normalizing the data.');
    % Eliminate extreme outliers, i.e., any point that exceedes 5 times the
    % inter quantile range, by bounding them to that value.
    disp('-> Removing extreme outliers from the feature values.');
    himask = bsxfun(@gt,X,model.bound.hibound);
    lomask = bsxfun(@lt,X,model.bound.lobound);
    X = X.*~(himask | lomask) + bsxfun(@times,himask,model.bound.hibound) + ...
                                bsxfun(@times,lomask,model.bound.lobound);
end

if model.opts.norm.flag
    % Normalize the data using Box-Cox and Z-transformations
    disp('-> Auto-normalizing the data.');
    X = bsxfun(@minus,X,model.norm.minX)+1;
    X = bsxfun(@rdivide,bsxfun(@power,X,model.norm.lambdaX)-1,model.norm.lambdaX);
    X = bsxfun(@rdivide,bsxfun(@minus,X,model.norm.muX),model.norm.sigmaX);
    
    Y(Y==0) = eps; % Assumes that Y is always positive and higher than 1e-16
    Y = bsxfun(@rdivide,bsxfun(@power,Y,model.norm.lambdaY)-1,model.norm.lambdaY);
    Y = bsxfun(@rdivide,bsxfun(@minus,Y,model.norm.muY),model.norm.sigmaY);
end
% ---------------------------------------------------------------------
% This is the final subset of features.
X = X(:,model.featsel.idx);
featlabels = strrep(varlabels(isfeat),'feature_','');
featlabels = featlabels(model.featsel.idx);
algolabels = strrep(varlabels(isalgo),'algo_','');
% ---------------------------------------------------------------------
%  Calculate the two dimensional projection using the PBLDR algorithm
%  (Munoz et al. Mach Learn 2018)
Z = X*model.pbldr.A';
% -------------------------------------------------------------------------
% Algorithm selection. Fit a model that would separate the space into
% classes of good and bad performance. 
Yhat = 0.*Ybin;
probs = 0.*Ybin;
for i=1:nalgos
    [Yhat(:,i),~,probs(:,i)] = svmpredict(randi([1 2], ninst,1), Z, model.algosel.svm{i}, '-q');
end
Yhat = Yhat==2; % Make it binary

end