function [X, out] = SIFTED2(X, Y, Ybin, featlabels, opts)
% SIFTED2  Automated feature selection for instance-space analysis.
%
%   [X, out] = SIFTED2(X, Y, Ybin, featlabels, opts)
%
%   Selects features in two stages:
%     1. Correlation filter  -- keeps features whose Pearson correlation
%        with any algorithm's performance exceeds opts.rho and is
%        statistically significant at the opts.pval level.
%     2. Correlation clustering + GA  -- groups remaining features by
%        correlation distance into opts.K clusters, then uses a genetic
%        algorithm with a KNN fitness function to pick one representative
%        feature per cluster.
%
%   Inputs
%     X          - (ninst x nfeats) feature matrix
%     Y          - (ninst x nalgos) performance matrix
%     Ybin       - (ninst x nalgos) logical good-performance matrix
%     featlabels - (1 x nfeats) cell array of feature name strings
%     opts       - struct with fields (see ISAdefaults for defaults):
%                    rho        double  minimum absolute correlation to keep (0.10)
%                    pval       double  significance level for correlation (0.05)
%                    K          int     number of feature clusters (10)
%                    MaxIter    int     k-means maximum iterations (1000)
%                    Replicates int     k-means replicates (100)
%
%   Outputs
%     X    - (ninst x nselected) reduced feature matrix
%     out  - struct with fields:
%              selvars   indices of selected features (into original X columns)
%              rho, p    full correlation and p-value matrices
%              eva       evalclusters result object
%              clust     (nfeats x K) cluster membership matrix
%              Ksuggested  (optional) better K value if current K gives poor silhouette
%
%   By: Mario Andres Munoz Acosta
%       School of Mathematics and Statistics
%       The University of Melbourne
%       Australia
%       2021

narginchk(5, 5);
if size(X,1) ~= size(Y,1) || size(X,1) ~= size(Ybin,1)
    error('ISA:SIFTED2:sizeMismatch', ...
        'X, Y and Ybin must have the same number of rows.');
end
if size(X,2) ~= length(featlabels)
    error('ISA:SIFTED2:labelMismatch', ...
        'featlabels length (%d) must match number of features in X (%d).', ...
        length(featlabels), size(X,2));
end
if ~isfield(opts, 'pval'), opts.pval = 0.05; end

% -------------------------------------------------------------------------
% GA and clustering parameters (fixed internal constants; not user-facing
% because they control optimiser behaviour, not the feature-selection logic)
Kfolds            = 5;    % CV folds for KNN fitness evaluation inside GA
FitnessLimit      = 0;    % GA stops when fitness <= this (0 = perfect KNN loss)
FunctionTolerance = 1e-3; % GA stops when improvement < this across stall window
MaxGenerations    = 100;  % hard cap on GA generations
MaxStallGenerations = 5;  % early-stop if best fitness unchanged for this many generations
PopulationSize    = 50;   % GA population size (caps combinations at ~50*K evaluations)
% Silhouette thresholds for advisory messages only
unnaceptableClustering = 0.50;
acceptableClustering   = 0.75;
% -------------------------------------------------------------------------

clearCache(); % reset the persistent fitness cache at the start of each call

if exist('gcp','file')==2
    mypool = gcp('nocreate');
    if ~isempty(mypool)
        nworkers = mypool.NumWorkers;
    else
        nworkers = 0;
    end
else
    nworkers = 0;
end

% -------------------------------------------------------------------------
nfeats = size(X, 2);
if nfeats <= 1
    error('ISA:SIFTED2:tooFewFeatures', ...
        'There is only 1 feature. Stopping space construction.');
elseif nfeats <= 3
    fprintf('[SIFTED] There are 3 or less features. Skipping feature selection.\n');
    out.selvars = 1:nfeats;
    return;
end
% -------------------------------------------------------------------------
fprintf('[SIFTED] Selecting features based on correlation with performance.\n');
[out.rho, out.p] = corr(X, Y, 'rows', 'pairwise');
rho = out.rho;
rho(isnan(rho) | (out.p > opts.pval)) = 0;
[rho, row] = sort(abs(rho), 1, 'descend');
out.selvars = false(1, nfeats);
% Always keep the most correlated feature for each algorithm
out.selvars(unique(row(1,:))) = true;
% Keep any feature with correlation at least opts.rho for any algorithm
for ii = 2:nfeats
    out.selvars(unique(row(ii, rho(ii,:) >= opts.rho))) = true;
end
out.selvars = find(out.selvars);
Xaux = X(:, out.selvars);
fprintf('[SIFTED] Keeping %d out of %d features (correlation).\n', size(Xaux,2), nfeats);
% -------------------------------------------------------------------------
nfeats = size(Xaux, 2);
if nfeats <= 1
    error('ISA:SIFTED2:tooFewFeatures', ...
        'There is only 1 feature after correlation filter. Stopping space construction.');
elseif nfeats <= 3
    fprintf('[SIFTED] There are 3 or less features. Skipping correlation clustering selection.\n');
    X = Xaux;
    return;
elseif nfeats <= opts.K
    fprintf('[SIFTED] Fewer features than clusters (%d <= %d). Skipping correlation clustering selection.\n', ...
        nfeats, opts.K);
    X = Xaux;
    return;
end
% -------------------------------------------------------------------------
fprintf('[SIFTED] Selecting features based on correlation clustering.\n');
state = rng;
rng('default');
out.eva = evalclusters(Xaux', 'kmeans', 'Silhouette', 'KList', 3:nfeats, ...
                              'Distance', 'correlation');
fprintf('[SIFTED] Average silhouette values for each number of clusters.\n');
disp([out.eva.InspectedK; out.eva.CriterionValues]);
if out.eva.CriterionValues(out.eva.InspectedK==opts.K) < unnaceptableClustering
    fprintf('[SIFTED] The silhouette value for K=%d is below %.2f. You should consider increasing K.\n', ...
        opts.K, unnaceptableClustering);
    out.Ksuggested = out.eva.InspectedK(find(out.eva.CriterionValues > acceptableClustering, 1));
    if ~isempty(out.Ksuggested)
        fprintf('[SIFTED] A suggested value of K is %d\n', out.Ksuggested);
    end
end
% -------------------------------------------------------------------------
rng('default');
out.clust = bsxfun(@eq, kmeans(Xaux', opts.K, 'Distance', 'correlation', ...
                                              'MaxIter', opts.MaxIter, ...
                                              'Replicates', opts.Replicates, ...
                                              'Options', statset('UseParallel', nworkers~=0), ...
                                              'OnlinePhase', 'on'), 1:opts.K);
rng(state);
fprintf('[SIFTED] Constructing %d clusters of features.\n', opts.K);
fprintf('[SIFTED] Using a GA+LookUpTable to find an optimal combination.\n');
% -------------------------------------------------------------------------
cvpart  = cvpartition(size(Xaux,1), 'Kfold', Kfolds);
fcnwrap = @(x) costfcn(x, Xaux, Y, Ybin, out.clust, cvpart, featlabels(out.selvars), nworkers);
gaopts  = optimoptions('ga', ...
    'FitnessLimit',      FitnessLimit, ...
    'FunctionTolerance', FunctionTolerance, ...
    'MaxGenerations',    MaxGenerations, ...
    'MaxStallGenerations', MaxStallGenerations, ...
    'PopulationSize',    PopulationSize);
ind = ga(fcnwrap, opts.K, [], [], [], [], ones(1,opts.K), sum(out.clust), [], 1:opts.K, gaopts);

decoder = false(1, size(Xaux,2));
for i = 1:opts.K
    aux = find(out.clust(:,i));
    decoder(aux(ind(i))) = true;
end
out.selvars = out.selvars(decoder);
X = X(:, out.selvars);
fprintf('[SIFTED] Keeping %d out of %d features (clustering).\n', size(X,2), nfeats);

end
% =========================================================================
function y = costfcn(ind, X, Y, Ybin, clust, cvpart, featlabels, nworkers)
    persistent mymap
    if isempty(mymap)
        mymap = containers.Map('KeyType','char','ValueType','double');
    end
    % Internal PILOT call uses minimal settings: 5 random starts, numerical
    % solver only. KNN uses 3 neighbours (D+1 for 2D projection).
    ntries      = 5;
    analytic    = false;
    kneighbours = 3;

    idx = false(1, size(X,2));
    for i = 1:length(ind)
        aux = find(clust(:,i));
        idx(aux(ind(i))) = true;
    end

    key = strrep(num2str(idx), ' ', '');
    if isKey(mymap, key)
        y = values(mymap, {key});
        y = y{1};
    else
        out = PILOT(X(:,idx), Y, featlabels(idx), struct('analytic', analytic, 'ntries', ntries));
        Z = out.Z;
        y = -Inf;
        % NOTE: ga runs fitness calls serially (no UseParallel in gaopts), so this
        % inner parfor does not nest and can safely use nworkers.
        parfor (ii = 1:size(Y,2), nworkers)
            knn = fitcknn(Z, Ybin(:,ii), 'CVPartition', cvpart, 'NumNeighbors', kneighbours);
            y = max(y, knn.kfoldLoss);
        end
        mymap(key) = y;
    end
end
% =========================================================================
function clearCache()
% Clears the persistent fitness cache in costfcn.
% Called at the start of each SIFTED2 invocation to prevent stale hits
% across successive buildIS calls in the same MATLAB session.
    clear costfcn
end
% =========================================================================
