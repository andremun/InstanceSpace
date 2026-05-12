function out = PYTHIA(Z, Y, Ybin, Ybest, algolabels, opts)
% -------------------------------------------------------------------------
% PYTHIA.m
% -------------------------------------------------------------------------
%
% By: Mario Andres Munoz Acosta
%     School of Mathematics and Statistics
%     The University of Melbourne
%     Australia
%     2020
%
% Unified binary classifier for instance-space algorithm selection.
% Classifier type is selected via opts.classifier ('knn' or 'svm').
% Hyperparameters are tuned with a Sobol quasi-random search + k-fold CV.
%
% -------------------------------------------------------------------------

% Handle deprecated flags silently for backward compatibility.
if isfield(opts, 'uselibsvm') && opts.uselibsvm
    warning('ISA:PYTHIA:libsvmDeprecated', ...
        ['opts.uselibsvm is deprecated and ignored. ' ...
         'Migrate saved models with ISAmigrateModel.']);
end
if isfield(opts, 'useknn') && opts.useknn && ~isfield(opts, 'classifier')
    opts.classifier = 'knn';
end

fprintf('  -> Initializing PYTHIA.\n');
[Znorm, out.mu, out.sigma] = zscore(Z);
[ninst, nalgos] = size(Ybin);
classifierType = opts.classifier;
nsobol         = opts.nsobol;

if opts.useweights
    fprintf('  -> Using cost-sensitive classification.\n');
    W = abs(Y - nanmean(Y(:)));
    W(W==0) = min(W(W~=0));
    W(isnan(W)) = max(W(~isnan(W)));
    out.W = W;
else
    fprintf('  -> Not using cost-sensitive classification.\n');
    W = ones(ninst, nalgos);
end

fprintf('-------------------------------------------------------------------------\n');
if strcmp(classifierType, 'svm')
    if opts.ispolykrnl
        kname = 'polynomial';
    else
        kname = 'gaussian';
    end
    fprintf('  -> Using a SVM classifier (%s kernel) with %d-point Sobol search.\n', kname, nsobol);
else
    fprintf('  -> Using a KNN classifier with %d-point Sobol search.\n', nsobol);
end
fprintf('  -> Using a %d-fold stratified cross-validation.\n', opts.cvfolds);
fprintf('-------------------------------------------------------------------------\n');
fprintf('  -> Training has started. PYTHIA may take a while to complete...\n');

out.cp             = cell(1, nalgos);
out.classifier     = cell(1, nalgos);
out.cvcmat         = zeros(nalgos, 4);
out.Ysub           = false(ninst, nalgos);
out.Yhat           = false(ninst, nalgos);
out.Pr0sub         = zeros(ninst, nalgos);
out.Pr0hat         = zeros(ninst, nalgos);
out.param1         = zeros(1, nalgos);
out.param2         = zeros(1, nalgos);
out.classifierType = classifierType;

t = tic;
for i = 1:nalgos
    tic;
    state = rng; rng('default');
    out.cp{i} = cvpartition(Ybin(:,i), 'Kfold', opts.cvfolds, 'Stratify', true);
    [out.classifier{i}, out.Ysub(:,i), out.Pr0sub(:,i), ...
     out.Yhat(:,i), out.Pr0hat(:,i), ...
     out.param1(i), out.param2(i)] = ...
        fitSobolClassifier(classifierType, Znorm, Ybin(:,i), W(:,i), ...
                           out.cp{i}, nsobol, opts.ispolykrnl);
    rng(state);
    out.cvcmat(i,:) = buildCvConfmat(Ybin(:,i), out.Ysub(:,i));
    if opts.verbose
        remaining = nalgos - i;
        if remaining == 0
            fprintf('    -> PYTHIA has trained a model for ''%s'', there are no models left to train.\n', algolabels{i});
        elseif remaining == 1
            fprintf('    -> PYTHIA has trained a model for ''%s'', there is 1 model left to train.\n', algolabels{i});
        else
            fprintf('    -> PYTHIA has trained a model for ''%s'', there are %d models left to train.\n', algolabels{i}, remaining);
        end
        fprintf('      -> Elapsed time: %.2fs\n', toc);
    end
end

tn = out.cvcmat(:,1); fp = out.cvcmat(:,3);
fn = out.cvcmat(:,2); tp = out.cvcmat(:,4);
out.precision = tp ./ (tp + fp);
out.recall    = tp ./ (tp + fn);
out.accuracy  = (tp + tn) ./ ninst;

fprintf('-------------------------------------------------------------------------\n');
fprintf('  -> PYTHIA has completed training the models.\n');
fprintf('  -> The average cross validated precision is: %.1f%%\n', 100.*nanmean(out.precision));
fprintf('  -> The average cross validated accuracy is: %.1f%%\n',  100.*nanmean(out.accuracy));
fprintf('      -> Elapsed time: %.2fs\n', toc(t));
fprintf('-------------------------------------------------------------------------\n');

% Algorithm selection: prefer the classifier with the best CV precision.
if nalgos > 1
    [best, out.selection0] = max(bsxfun(@times, out.Yhat, out.precision'), [], 2);
else
    best = out.Yhat;
    out.selection0 = double(out.Yhat);
end
[~, default] = max(mean(Ybin));
out.selection1 = out.selection0;
out.selection0(best <= 0) = 0;
out.selection1(best <= 0) = default;

sel0 = bsxfun(@eq, out.selection0, 1:nalgos);
sel1 = bsxfun(@eq, out.selection1, 1:nalgos);
avgperf = nanmean(Y);
stdperf = nanstd(Y);
Yfull = Y; Ysvms = Y;
Y(~sel0) = NaN;
Yfull(~sel1) = NaN;
Ysvms(~out.Yhat) = NaN;

pgood = mean(any(Ybin & sel1, 2));
fb = sum(any( Ybin & ~sel0, 2));
fg = sum(any(~Ybin &  sel0, 2));
tg = sum(any( Ybin &  sel0, 2));
precisionsel = tg ./ (tg + fg);
recallsel    = tg ./ (tg + fb);

% Summary column labels for hyperparameters depend on classifier type.
if strcmp(classifierType, 'knn')
    p1label = 'NumNeighbors';
    p2label = 'Distance';
    distOpts = {'euclidean','cityblock','cosine','correlation'};
    p2cells = cellfun(@(x) distOpts{min(4,max(1,round(x)))}, ...
                      num2cell(out.param2), 'UniformOutput', false);
else
    p1label = 'BoxConstraint';
    p2label = 'KernelScale';
    p2cells = num2cell(round(out.param2, 3));
end

fprintf('  -> PYTHIA is preparing the summary table.\n');
out.summary = cell(nalgos+3, 11);
out.summary{1,1} = 'Algorithms ';
out.summary(2:end-2, 1) = algolabels;
out.summary(end-1:end, 1) = {'Oracle','Selector'};
out.summary(1, 2:11) = {'Avg_Perf_all_instances';
                         'Std_Perf_all_instances';
                         'Probability_of_good';
                         'Avg_Perf_selected_instances';
                         'Std_Perf_selected_instances';
                         'CV_model_accuracy';
                         'CV_model_precision';
                         'CV_model_recall';
                         p1label;
                         p2label};
out.summary(2:end, 2)    = num2cell(round([avgperf nanmean(Ybest) nanmean(Yfull(:))], 3));
out.summary(2:end, 3)    = num2cell(round([stdperf nanstd(Ybest) nanstd(Yfull(:))], 3));
out.summary(2:end, 4)    = num2cell(round([mean(Ybin) 1 pgood], 3));
out.summary(2:end, 5)    = num2cell(round([nanmean(Ysvms) NaN nanmean(Y(:))], 3));
out.summary(2:end, 6)    = num2cell(round([nanstd(Ysvms) NaN nanstd(Y(:))], 3));
out.summary(2:end, 7)    = num2cell(round(100.*[out.accuracy' NaN NaN], 1));
out.summary(2:end, 8)    = num2cell(round(100.*[out.precision' NaN precisionsel], 1));
out.summary(2:end, 9)    = num2cell(round(100.*[out.recall' NaN recallsel], 1));
out.summary(2:end-2, 10) = num2cell(round(out.param1, 3));
out.summary(2:end-2, 11) = p2cells;
out.summary(cellfun(@(x) isnumeric(x) && all(isnan(x(:))), out.summary)) = {[]};
fprintf('  -> PYTHIA has completed! Performance of the models:\n\n');
disp(out.summary);
end

% =========================================================================
% SUBFUNCTIONS
% =========================================================================

function [clf, Ysub, Psub, Yhat, Phat, p1, p2] = ...
        fitSobolClassifier(type, Z, Ybin, W, cp, nsobol, ispolykrnl)
% Tune hyperparameters via Sobol search with k-fold CV, then train final model.
ninst    = size(Z, 1);
nworkers = getParallelWorkers();

% Generate nsobol quasi-random candidates in [0,1]^2.
ss = sobolset(2, 'Skip', 1);
X  = net(ss, nsobol);

params = zeros(nsobol, 2);
if strcmp(type, 'svm')
    logRng = [-10, 4];
    params(:,1) = 2.^(logRng(1) + X(:,1) .* diff(logRng));  % BoxConstraint
    params(:,2) = 2.^(logRng(1) + X(:,2) .* diff(logRng));  % KernelScale
else  % knn
    params(:,1) = max(1, round(1 + X(:,1) .* 24));           % NumNeighbors [1..25]
    params(:,2) = max(1, min(4, ceil(X(:,2) .* 4)));          % Distance index [1..4]
end

% Evaluate all Sobol candidates via k-fold CV.
Ysub_all = false(ninst, nsobol);
Psub_all = zeros(ninst, nsobol);

for fold = 1:cp.NumTestSets
    itrain = logical(cp.training(fold));
    itest  = logical(cp.test(fold));
    Ztrain = Z(itrain, :);
    Ytrain = logical(Ybin(itrain));
    Wtrain = W(itrain);
    Ztest  = Z(itest, :);
    ntest  = sum(itest);

    Yfold = false(ntest, nsobol);
    Pfold = zeros(ntest, nsobol);
    parfor (j = 1:nsobol, nworkers)
        [Yfold(:,j), Pfold(:,j)] = evalFoldClassifier( ...
            type, Ztrain, Ytrain, Wtrain, Ztest, params(j,:), ispolykrnl);
    end
    Ysub_all(itest, :) = Yfold;
    Psub_all(itest, :) = Pfold;
end

% Select candidate with minimum CV misclassification error.
Ybin_rep = repmat(logical(Ybin), 1, nsobol);
errs = mean(Ysub_all ~= Ybin_rep, 1);
[~, best] = min(errs);
Ysub = Ysub_all(:, best);
Psub = Psub_all(:, best);
p1   = params(best, 1);
p2   = params(best, 2);

[clf, Yhat, Phat] = trainFinalClassifier( ...
    type, Z, logical(Ybin), W, params(best,:), ispolykrnl);
end

% -------------------------------------------------------------------------

function [Yp, Pp] = evalFoldClassifier(type, Ztrain, Ytrain, Wtrain, Ztest, params_j, ispolykrnl)
% Train one classifier on a CV fold and return predictions on the test fold.
distOpts = {'euclidean','cityblock','cosine','correlation'};
if strcmp(type, 'svm')
    if ispolykrnl; kernel = 'polynomial'; else; kernel = 'gaussian'; end
    clf = fitcsvm(Ztrain, Ytrain, ...
                  'Weights', Wtrain, ...
                  'BoxConstraint', params_j(1), ...
                  'KernelScale',   params_j(2), ...
                  'KernelFunction', kernel, ...
                  'Standardize', false, ...
                  'CacheSize', 'maximal');
    try
        clf = fitSVMPosterior(clf);
        [Yp, aux] = clf.predict(Ztest);
        Pp = aux(:,1);
    catch
        [Yp, ~] = clf.predict(Ztest);
        Pp = zeros(size(Ztest,1), 1);
    end
else  % knn
    dist = distOpts{params_j(2)};
    clf = fitcknn(Ztrain, Ytrain, ...
                  'Weights', Wtrain, ...
                  'NumNeighbors', round(params_j(1)), ...
                  'Distance', dist, ...
                  'Standardize', false);
    [Yp, aux] = clf.predict(Ztest);
    Pp = aux(:,1);
end
Yp = logical(Yp);
end

% -------------------------------------------------------------------------

function [clf, Yhat, Phat] = trainFinalClassifier(type, Z, Ybin, W, params_best, ispolykrnl)
% Train the final model on all data with the best hyperparameters.
distOpts = {'euclidean','cityblock','cosine','correlation'};
if strcmp(type, 'svm')
    if ispolykrnl; kernel = 'polynomial'; else; kernel = 'gaussian'; end
    clf = fitcsvm(Z, Ybin, ...
                  'Weights', W, ...
                  'BoxConstraint', params_best(1), ...
                  'KernelScale',   params_best(2), ...
                  'KernelFunction', kernel, ...
                  'Standardize', false, ...
                  'CacheSize', 'maximal', ...
                  'RemoveDuplicates', true);
    try
        clf = fitSVMPosterior(clf);
        [Yhat, aux] = clf.predict(Z);
        Phat = aux(:,1);
    catch
        [Yhat, ~] = clf.predict(Z);
        Phat = zeros(size(Z,1), 1);
    end
else  % knn
    dist = distOpts{params_best(2)};
    clf = fitcknn(Z, Ybin, ...
                  'Weights', W, ...
                  'NumNeighbors', round(params_best(1)), ...
                  'Distance', dist, ...
                  'Standardize', false);
    [Yhat, aux] = clf.predict(Z);
    Phat = aux(:,1);
end
Yhat = logical(Yhat);
end

% -------------------------------------------------------------------------

function cm = buildCvConfmat(Ytrue, Ypred)
% Build a [TN FN FP TP] row vector from binary labels (handles degenerate 1x1 case).
aux = confusionmat(logical(Ytrue), logical(Ypred));
if numel(aux) ~= 4
    caux = aux(1);
    aux  = zeros(2);
    if all(Ytrue == 0)
        if all(Ypred == 0)
            aux(1,1) = caux;
        elseif all(Ypred == 1)
            aux(2,1) = caux;
        end
    elseif all(Ytrue == 1)
        if all(Ypred == 0)
            aux(1,2) = caux;
        elseif all(Ypred == 1)
            aux(2,2) = caux;
        end
    end
end
cm = aux(:)';
end

% -------------------------------------------------------------------------

function nw = getParallelWorkers()
nw = 0;
if exist('gcp', 'file') == 2
    pool = gcp('nocreate');
    if ~isempty(pool)
        nw = pool.NumWorkers;
    end
end
end
% =========================================================================
