function out = PYTHIA(Z, Y, Ybin, Ybest, algolabels, opts, trainedModel)
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
%
% Training mode  (6 args):
%   out = PYTHIA(Z, Y, Ybin, Ybest, algolabels, opts)
%
% Evaluation mode (7 args — uses classifiers trained by a previous buildIS):
%   out = PYTHIA(Z, Y, Ybin, Ybest, algolabels, opts, trainedModel)
%
% Classifier is selected via opts.classifier (default 'knn').
% See ISAgetClassifierFcn for the full registry.
%
% Hyperparameter tuning (opts.tuning):
%   'sobol' (default) — scrambled Sobol quasi-random search, opts.nTuningIter evals
%   'none'            — use pre-supplied opts.params directly; skip tuning
%   'bayes'           — reserved; falls back to 'sobol' with a warning
% -------------------------------------------------------------------------

% Defensive defaults — guard against callers that bypass ISAdefaults.
if ~isfield(opts, 'classifier'),      opts.classifier      = 'knn';    end
if ~isfield(opts, 'tuning'),          opts.tuning          = 'sobol';  end
if ~isfield(opts, 'nTuningIter'),     opts.nTuningIter     = 20;       end
if ~isfield(opts, 'cvfolds'),         opts.cvfolds         = 5;        end
if ~isfield(opts, 'ispolykrnl'),      opts.ispolykrnl      = false;    end
if ~isfield(opts, 'useweights'),      opts.useweights      = false;    end
if ~isfield(opts, 'verbose'),         opts.verbose         = true;     end
if ~isfield(opts, 'params'),          opts.params          = [];       end
if ~isfield(opts, 'skip'),            opts.skip            = false;    end
if ~isfield(opts, 'ensembleMethod'),  opts.ensembleMethod  = 'Bag';    end
if ~isfield(opts, 'seed'),            opts.seed            = 42;       end
% Handle deprecated flags.
if isfield(opts, 'uselibsvm') && opts.uselibsvm
    warning('ISA:PYTHIA:libsvmDeprecated', ...
        ['opts.uselibsvm is deprecated and ignored. ' ...
         'Migrate saved models with ISAmigrateModel.']);
end
if strcmp(opts.tuning, 'bayes')
    warning('ISA:PYTHIA:bayesFallback', ...
        ['opts.tuning=''bayes'' is not yet fully supported; ' ...
         'falling back to ''sobol''.']);
    opts.tuning = 'sobol';
end

% -------------------------------------------------------------------------
% Eval mode — apply trained classifiers to new data, skip training.
if nargin == 7
    out = PYTHIAevalMode(Z, Y, Ybin, Ybest, algolabels, opts, trainedModel);
    return;
end

% -------------------------------------------------------------------------
% Training mode
% -------------------------------------------------------------------------
% Validate classifier name early (gives clean error before any work).
[~, p1label, p2label] = ISAgetClassifierFcn(opts.classifier);

fprintf('  -> Initializing PYTHIA.\n');
[Znorm, out.mu, out.sigma] = zscore(Z);
[ninst, nalgos] = size(Ybin);
classifierType = opts.classifier;

% Skip mode: fill outputs with NaN and return so TRACE can use Ybin directly.
if opts.skip
    fprintf('  -> opts.skip is true; bypassing classifier training.\n');
    mu_z = out.mu; sigma_z = out.sigma;      % preserve correct feature-wise params
    out = emptyPYTHIAout(ninst, nalgos, algolabels, Y, Ybin, Ybest, size(Z,2));
    out.mu = mu_z; out.sigma = sigma_z;
    return;
end

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

nIter = opts.nTuningIter;
precalcparams = ~isempty(opts.params) && isnumeric(opts.params) && ...
                size(opts.params,1)==nalgos && size(opts.params,2)==2;
if strcmp(opts.tuning, 'none') && ~precalcparams
    error('ISA:PYTHIA:noParamsForNoneTuning', ...
        ['opts.tuning=''none'' requires opts.params to be a valid [%d x 2] numeric ' ...
         'matrix of pre-calculated hyperparameters. ' ...
         'Either supply opts.params or change opts.tuning to ''sobol''.'], nalgos);
end

fprintf('-------------------------------------------------------------------------\n');
fprintf('  -> Classifier: %s | Tuning: %s (%d evals) | CV: %d-fold.\n', ...
        classifierType, opts.tuning, nIter, opts.cvfolds);
fprintf('-------------------------------------------------------------------------\n');
fprintf('  -> Training has started. PYTHIA may take a while to complete...\n');

out.classifiers    = cell(1, nalgos);
out.classifierType = classifierType;
out.cp             = cell(1, nalgos);
out.cvcmat         = zeros(nalgos, 4);
out.Ysub           = false(ninst, nalgos);
out.Yhat           = false(ninst, nalgos);
out.Pr0sub         = zeros(ninst, nalgos);
out.Pr0hat         = zeros(ninst, nalgos);
out.param1         = zeros(1, nalgos);
out.param2         = zeros(1, nalgos);
out.param2Label    = cell(1, nalgos);  % human-readable param2 (distance name for KNN)

% Save global RNG state so downstream code (TRACE, user scripts) is unaffected.
% Per-algorithm seeding inside the loop keeps each classifier's training
% reproducible without permanently mutating the session RNG.
prevRNG = rng;
t = tic;
for i = 1:nalgos
    tic;
    % Per-algorithm reproducible seed: opts.seed + i (per spec §6.2).
    rng(opts.seed + i, 'twister');
    out.cp{i} = cvpartition(Ybin(:,i), 'Kfold', opts.cvfolds, 'Stratify', true);

    if precalcparams
        p1_best = opts.params(i,1);
        p2_best = opts.params(i,2);
        % CV predictions with the pre-supplied params.
        [out.Ysub(:,i), out.Pr0sub(:,i)] = crossValPredict( ...
            classifierType, Znorm, Ybin(:,i), W(:,i), out.cp{i}, ...
            p1_best, p2_best, opts);
    else
        % Scrambled Sobol search.
        ss = sobolset(2, 'Skip', 1, 'Scramble', 'MatousekAffineOwen');
        X  = net(ss, nIter);
        [P1, P2] = sobolToParams(classifierType, X);

        [out.Ysub(:,i), out.Pr0sub(:,i), p1_best, p2_best] = ...
            sobolSearch(classifierType, Znorm, Ybin(:,i), W(:,i), ...
                        out.cp{i}, P1, P2, opts);
    end

    [out.classifiers{i}, out.Yhat(:,i), out.Pr0hat(:,i)] = ...
        trainFinalClassifier(classifierType, Znorm, logical(Ybin(:,i)), W(:,i), ...
                             p1_best, p2_best, opts);
    out.param1(i) = p1_best;
    out.param2(i) = p2_best;
    if strcmpi(classifierType, 'knn')
        distOpts = {'euclidean','cityblock','cosine','correlation'};
        out.param2Label{i} = distOpts{min(4,max(1,round(p2_best)))};
    end
    % Non-KNN classifiers: param2 is already numeric; param2Label left empty.

    cm = confusionmat(logical(Ybin(:,i)), out.Ysub(:,i), 'Order', [false true]);
    out.cvcmat(i,:) = cm(:)';

    if opts.verbose
        remaining = nalgos - i;
        if remaining == 0
            fprintf('    -> PYTHIA trained ''%s''; no models left.\n', algolabels{i});
        elseif remaining == 1
            fprintf('    -> PYTHIA trained ''%s''; 1 model left.\n', algolabels{i});
        else
            fprintf('    -> PYTHIA trained ''%s''; %d models left.\n', algolabels{i}, remaining);
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
fprintf('  -> PYTHIA training complete.\n');
fprintf('  -> Average CV precision: %.1f%%\n', 100.*nanmean(out.precision));
fprintf('  -> Average CV accuracy : %.1f%%\n', 100.*nanmean(out.accuracy));
fprintf('      -> Elapsed time: %.2fs\n', toc(t));
fprintf('-------------------------------------------------------------------------\n');

rng(prevRNG);  % restore global RNG so downstream code is unaffected

out = computeSelection(out, nalgos, Ybin);
out.summary = buildSummary(out, algolabels, nalgos, ninst, ...
                           Y, Ybin, Ybest, p1label, p2label);
fprintf('  -> PYTHIA has completed! Performance of the models:\n\n');
disp(out.summary);
end

% =========================================================================
%  EVAL MODE
% =========================================================================
function out = PYTHIAevalMode(Z, Y, Ybin, Ybest, algolabels, opts, trained)
% Apply trained classifiers to new (test) data.

if isfield(trained, 'mu') && isfield(trained, 'sigma')
    Znorm = (Z - trained.mu) ./ trained.sigma;
else
    Znorm = Z;
end

% Dispatch on field: new 'classifiers' (plural) or legacy 'svm' / 'knn'.
if isfield(trained, 'classifiers')
    clfs = trained.classifiers;
elseif isfield(trained, 'svm')
    clfs = trained.svm;
elseif isfield(trained, 'knn')
    clfs = trained.knn;
else
    error('ISA:PYTHIA:noClassifier', ...
        'Trained model has no recognised classifier field (classifiers / svm / knn).');
end

nalgos = length(clfs);
ninst  = size(Znorm, 1);
Y      = Y(:,1:nalgos);
Ybin   = Ybin(:,1:nalgos);

out.Yhat   = false(ninst, nalgos);
out.Pr0hat = zeros(ninst, nalgos);
out.cvcmat = zeros(nalgos, 4);

for ii = 1:nalgos
    clf = clfs{ii};
    if isempty(clf)
        continue;
    end
    if isstruct(clf)
        % Legacy LIBSVM struct.
        Yin = double(Ybin(:,ii)) + 1;
        [aux, ~, out.Pr0hat(:,ii)] = svmpredict(Yin, Znorm, clf, '-q');
        out.Yhat(:,ii) = logical(aux == 2);
    else
        [out.Yhat(:,ii), aux] = clf.predict(Znorm);
        out.Yhat(:,ii) = logical(out.Yhat(:,ii));
        if size(aux,2) >= 1; out.Pr0hat(:,ii) = aux(:,1); end
    end
    cm = confusionmat(logical(Ybin(:,ii)), out.Yhat(:,ii), 'Order', [false true]);
    out.cvcmat(ii,:) = cm(:)';
end

tn = out.cvcmat(:,1); fp = out.cvcmat(:,3);
fn = out.cvcmat(:,2); tp = out.cvcmat(:,4);
out.precision = tp ./ (tp + fp);
out.recall    = tp ./ (tp + fn);
out.accuracy  = (tp + tn) ./ ninst;

% Use the training-time precision for selection weighting if available; fall back
% to the freshly computed eval precision so old/migrated models still work.
if isfield(trained, 'precision') && ~isempty(trained.precision)
    out = computeSelection(out, nalgos, Ybin, trained.precision);
else
    out = computeSelection(out, nalgos, Ybin);
end
% Eval-mode summary has 9 columns (no hyperparameter columns).
out.summary = buildSummary(out, algolabels, nalgos, ninst, ...
                           Y, Ybin, Ybest, [], []);
fprintf('  -> PYTHIA has completed! Performance of the models:\n\n');
disp(out.summary);
end

% =========================================================================
%  SUBFUNCTIONS
% =========================================================================

function [Ysub, Psub, p1_best, p2_best] = sobolSearch( ...
        type, Z, Ybin, W, cp, P1, P2, opts)
% Evaluate nIter Sobol candidates via k-fold CV and return best params + CV predictions.
nsobol = length(P1);
ninst  = size(Z, 1);
nworkers = getParallelWorkers();

Ysub_all = false(ninst, nsobol);
Psub_all = zeros(ninst, nsobol);

for fold = 1:cp.NumTestSets
    itrain = logical(cp.training(fold));
    itest  = logical(cp.test(fold));
    Ztrain = Z(itrain,:);  Ytrain = logical(Ybin(itrain));  Wtrain = W(itrain);
    Ztest  = Z(itest,:);   ntest  = sum(itest);

    Yfold = false(ntest, nsobol);
    Pfold = zeros(ntest, nsobol);
    parfor (j = 1:nsobol, nworkers)
        [Yfold(:,j), Pfold(:,j)] = evalFoldClassifier( ...
            type, Ztrain, Ytrain, Wtrain, Ztest, P1(j), P2(j), opts);
    end
    Ysub_all(itest,:) = Yfold;
    Psub_all(itest,:) = Pfold;
end

Ybin_rep = repmat(logical(Ybin), 1, nsobol);
errs = mean(Ysub_all ~= Ybin_rep, 1);
% Invalidate candidates where any fold's training failed (NaN probability).
errs(any(isnan(Psub_all), 1)) = Inf;
if all(isinf(errs))
    warning('ISA:PYTHIA:allSobolFailed', ...
        'All Sobol candidates failed (training error on every fold). Using first candidate.');
    errs(1) = 0;
end
[~, best] = min(errs);
Ysub   = Ysub_all(:, best);
Psub   = Psub_all(:, best);
p1_best = P1(best);
p2_best = P2(best);
end

% -------------------------------------------------------------------------
function [Ysub, Psub] = crossValPredict(type, Z, Ybin, W, cp, p1, p2, opts)
% Run k-fold CV with fixed hyperparameters; return fold-level predictions.
ninst    = size(Z, 1);
nworkers = getParallelWorkers();
Ysub = false(ninst, 1);
Psub = zeros(ninst, 1);

for fold = 1:cp.NumTestSets
    itrain = logical(cp.training(fold));
    itest  = logical(cp.test(fold));
    Ztrain = Z(itrain,:);  Ytrain = logical(Ybin(itrain));  Wtrain = W(itrain);
    Ztest  = Z(itest,:);
    [Yfold, Pfold] = evalFoldClassifier(type, Ztrain, Ytrain, Wtrain, Ztest, p1, p2, opts);
    Ysub(itest) = Yfold;
    Psub(itest) = Pfold;
end
end

% -------------------------------------------------------------------------
function [Yp, Pp] = evalFoldClassifier(type, Ztrain, Ytrain, Wtrain, Ztest, p1, p2, opts)
% Train one classifier on a CV fold and predict on the test fold.
try
    clf = fitOneClassifier(type, Ztrain, Ytrain, Wtrain, p1, p2, opts);
    [Yp, Pp] = predictClassifier(clf, Ztest);
catch ME
    if opts.verbose
        warning('ISA:PYTHIA:foldTrainFailed', ...
            'Classifier training failed on CV fold (p1=%.4g, p2=%.4g): %s', ...
            p1, p2, ME.message);
    end
    Yp = false(size(Ztest,1), 1);
    Pp = NaN(size(Ztest,1), 1);  % NaN disqualifies this candidate from Sobol selection
end
end

% -------------------------------------------------------------------------
function [clf, Yhat, Phat] = trainFinalClassifier(type, Z, Ybin, W, p1, p2, opts)
% Train the final model on all data with the best hyperparameters.
clf = fitOneClassifier(type, Z, Ybin, W, p1, p2, opts, true);
[Yhat, Phat] = predictClassifier(clf, Z);
end

% -------------------------------------------------------------------------
function clf = fitOneClassifier(type, Z, Y, W, p1, p2, opts, isFinal)
% Dispatch training call to the appropriate MATLAB fitc* function.
if nargin < 8; isFinal = false; end

switch lower(type)
    case 'knn'
        distOpts = {'euclidean','cityblock','cosine','correlation'};
        clf = fitcknn(Z, Y, ...
                      'Weights', W, ...
                      'NumNeighbors', max(1, round(p1)), ...
                      'Distance', distOpts{max(1,min(4,round(p2)))}, ...
                      'Standardize', false);

    case 'svm'
        if opts.ispolykrnl; kfn = 'polynomial'; else; kfn = 'gaussian'; end
        args = {'Weights', W, ...
                'BoxConstraint', p1, 'KernelScale', p2, ...
                'KernelFunction', kfn, ...
                'Standardize', false, 'CacheSize', 'maximal'};
        if isFinal; args = [args, {'RemoveDuplicates', true}]; end
        clf = fitcsvm(Z, Y, args{:});
        try
            clf = fitSVMPosterior(clf);
        catch ME
            warning('ISA:PYTHIA:posteriorFailed', ...
                'fitSVMPosterior failed; predict will return raw decision scores: %s', ...
                ME.message);
        end

    case 'tree'
        clf = fitctree(Z, Y, ...
                       'Weights', W, ...
                       'MinLeafSize', max(1, round(p1)));

    case 'nb'
        try
            clf = fitcnb(Z, Y, 'Weights', W, ...
                         'DistributionNames', 'kernel', ...
                         'Width', p1);
        catch
            clf = fitcnb(Z, Y, 'Weights', W);
        end

    case 'linear'
        try
            clf = fitclinear(Z, Y, 'Weights', W, ...
                             'Lambda', p1, 'Learner', 'logistic');
        catch
            clf = fitclinear(Z, Y, 'Weights', W, 'Lambda', p1);
        end

    case 'ensemble'
        clf = fitcensemble(Z, Y, ...
                           'Weights', W, ...
                           'Method', opts.ensembleMethod, ...
                           'NumLearningCycles', max(10, round(p1)), ...
                           'Learners', templateTree('MinLeafSize', max(1, round(p2))));
end
end

% -------------------------------------------------------------------------
function [Yp, Pp] = predictClassifier(clf, Z)
[Yp, aux] = clf.predict(Z);
Yp = logical(Yp);
if size(aux,2) >= 1
    Pp = aux(:,1);
else
    Pp = zeros(size(Z,1), 1);
end
end

% -------------------------------------------------------------------------
function [P1, P2] = sobolToParams(type, X)
% Map Sobol points in [0,1]^2 to classifier-specific hyperparameter values.
switch lower(type)
    case 'knn'
        P1 = max(1, round(1 + X(:,1)*24));           % NumNeighbors: 1–25
        P2 = max(1, min(4, ceil(X(:,2)*4)));          % Distance index: 1–4
    case 'svm'
        P1 = 2.^(-10 + X(:,1)*14);                   % BoxConstraint: 2^[-10,4]
        P2 = 2.^(-10 + X(:,2)*14);                   % KernelScale: 2^[-10,4]
    case 'tree'
        P1 = max(1, round(1 + X(:,1)*99));            % MinLeafSize: 1–100
        P2 = ones(size(X,1), 1);                      % unused
    case 'nb'
        P1 = 10.^(-3 + X(:,1)*4);                    % Bandwidth: 1e-3–10
        P2 = ones(size(X,1), 1);                      % unused
    case 'linear'
        P1 = 10.^(-6 + X(:,1)*9);                    % Lambda: 1e-6–1e3
        P2 = ones(size(X,1), 1);                      % unused
    case 'ensemble'
        P1 = max(10, round(10 + X(:,1)*190));         % NumLearningCycles: 10–200
        P2 = max(1, round(1 + X(:,2)*19));            % MinLeafSize: 1–20
end
end

% -------------------------------------------------------------------------
function out = computeSelection(out, nalgos, Ybin, trainPrecision)
% Compute algorithm selection vectors using CV-precision-weighted voting.
if nargin < 4; trainPrecision = out.precision; end
% Replace NaN precision (tp+fp==0) with 0 so max() never returns NaN.
trainPrecision(isnan(trainPrecision)) = 0;
if nalgos > 1
    [best, out.selection0] = max(bsxfun(@times, out.Yhat, trainPrecision'), [], 2);
else
    best = out.Yhat;
    out.selection0 = double(out.Yhat);
end
[~, default] = max(mean(Ybin));
out.selection1 = out.selection0;
out.selection0(best <= 0) = 0;
out.selection1(best <= 0) = default;
end

% -------------------------------------------------------------------------
function summary = buildSummary(out, algolabels, nalgos, ninst, ...
                                Y, Ybin, Ybest, p1label, p2label)
% Build the PYTHIA summary cell array.
hasP2 = ~isempty(p1label) && ~strcmp(p2label, 'N/A');
ncols = 9 + ~isempty(p1label) + hasP2;  % 9 (eval), 10 (1-param classifier), 11 (2-param)

sel0 = bsxfun(@eq, out.selection0, 1:nalgos);
sel1 = bsxfun(@eq, out.selection1, 1:nalgos);
avgperf = nanmean(Y);
stdperf = nanstd(Y);
Yfull = Y; Ysvms = Y;
Y(~sel0)     = NaN;
Yfull(~sel1) = NaN;
Ysvms(~out.Yhat) = NaN;

pgood = mean(any(Ybin & sel1, 2));
fb = sum(any( Ybin & ~sel0, 2));
fg = sum(any(~Ybin &  sel0, 2));
tg = sum(any( Ybin &  sel0, 2));
precisionsel = tg ./ (tg + fg);
recallsel    = tg ./ (tg + fb);

summary = cell(nalgos+3, ncols);
summary{1,1} = 'Algorithms ';
summary(2:end-2, 1) = algolabels(1:nalgos);
summary(end-1:end, 1) = {'Oracle','Selector'};
colheads = {'Avg_Perf_all_instances';
            'Std_Perf_all_instances';
            'Probability_of_good';
            'Avg_Perf_selected_instances';
            'Std_Perf_selected_instances';
            'CV_model_accuracy';
            'CV_model_precision';
            'CV_model_recall'};
if ~isempty(p1label)
    colheads = [colheads; {p1label}];
    if hasP2; colheads = [colheads; {p2label}]; end
end
summary(1, 2:end) = colheads;
% Column vectors (nalgos+2)×1 assembled with ; — required for 2D cell-slice assignment.
summary(2:end, 2) = num2cell(round([avgperf(:);        nanmean(Ybest);  nanmean(Yfull(:))], 3));
summary(2:end, 3) = num2cell(round([stdperf(:);        nanstd(Ybest);   nanstd(Yfull(:))],  3));
summary(2:end, 4) = num2cell(round([mean(Ybin)';       1;               pgood],             3));
summary(2:end, 5) = num2cell(round([nanmean(Ysvms)';   NaN;             nanmean(Y(:))],     3));
summary(2:end, 6) = num2cell(round([nanstd(Ysvms)';    NaN;             nanstd(Y(:))],      3));
summary(2:end, 7) = num2cell(round(100.*[out.accuracy;  NaN;            NaN],               1));
summary(2:end, 8) = num2cell(round(100.*[out.precision; NaN;            precisionsel],       1));
summary(2:end, 9) = num2cell(round(100.*[out.recall;    NaN;            recallsel],          1));
if ~isempty(p1label)
    summary(2:end-2, 10) = num2cell(round(out.param1(:), 3));
    if hasP2
        % KNN: param2Label holds the resolved distance name (categorical string).
        % All other classifiers: use the numeric param2 value directly so the
        % summary/CSV cell contains a number, not a string representation.
        if strcmpi(out.classifierType, 'knn') && isfield(out, 'param2Label') ...
                && ~isempty(out.param2Label{1})
            summary(2:end-2, 11) = out.param2Label(:);
        else
            summary(2:end-2, 11) = num2cell(round(out.param2(:), 3));
        end
    end
end
summary(cellfun(@(x) isnumeric(x) && all(isnan(x(:))), summary)) = {[]};
end

% -------------------------------------------------------------------------
function out = emptyPYTHIAout(ninst, nalgos, algolabels, Y, Ybin, Ybest, nfeats)
% Return a no-op PYTHIA output when opts.skip = true.
% nfeats must match size(Z,2) so that eval-mode normalisation (Z - mu)./sigma works.
% mu/sigma are overwritten by the caller with the correct zscore parameters.
%
% precision/recall/accuracy are nalgos×1 column vectors (matching training/eval
% mode) so that buildSummary's transpose-and-concatenate idiom works correctly.
if nargin < 7; nfeats = 2; end   % 2D projected space is the common case
out.classifiers    = cell(1, nalgos);
out.classifierType = 'none';
out.mu             = zeros(1, nfeats);
out.sigma          = ones(1, nfeats);
out.cp             = cell(1, nalgos);
out.cvcmat         = zeros(nalgos, 4);
out.Ysub           = false(ninst, nalgos);
out.Yhat           = false(ninst, nalgos);
out.Pr0sub         = zeros(ninst, nalgos);
out.Pr0hat         = zeros(ninst, nalgos);
out.precision      = NaN(nalgos, 1);   % column vector — matches training/eval orientation
out.recall         = NaN(nalgos, 1);
out.accuracy       = NaN(nalgos, 1);
out.param1         = zeros(1, nalgos);
out.param2         = zeros(1, nalgos);
out.param2Label    = cell(1, nalgos);
out.selection0     = zeros(ninst, 1);
out.selection1     = zeros(ninst, 1);
% Populate a NaN-filled 9-column summary so scriptcsv/scriptpng can index
% into it without error even in skip mode.  No hyperparameter columns since
% no classifier was trained.
out.summary = buildSummary(out, algolabels, nalgos, ninst, Y, Ybin, Ybest, [], []);
end

% -------------------------------------------------------------------------
function nw = getParallelWorkers()
nw = 0;
if exist('gcp', 'file') == 2
    pool = gcp('nocreate');
    if ~isempty(pool); nw = pool.NumWorkers; end
end
end
% =========================================================================
