function [X, out] = SIFTED(X, Y, Ybin, featlabels, opts)
% SIFTED  Automated feature selection for instance-space analysis.
%
%   [X, out] = SIFTED(X, Y, Ybin, featlabels, opts)
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
%
%   Simpson, C., Munoz, M.A., Kandanaarachchi, S. & Campello, R.J.G.B.
%   (2025). ISA3: A 3-dimensional expansion of Instance Space Analysis.
%   Machine Learning, 114, 240. https://doi.org/10.1007/s10994-025-06871-5
% -------------------------------------------------------------------------

narginchk(5, 5);
if size(X,1) ~= size(Y,1) || size(X,1) ~= size(Ybin,1)
    error('ISA:SIFTED:sizeMismatch', ...
        'X, Y and Ybin must have the same number of rows.');
end
if size(X,2) ~= length(featlabels)
    error('ISA:SIFTED:labelMismatch', ...
        'featlabels length (%d) must match number of features in X (%d).', ...
        length(featlabels), size(X,2));
end
if ~isfield(opts, 'pval'), opts.pval = 0.05; end
if ~isfield(opts, 'dims'), opts.dims = 2;    end
if ~(isnumeric(opts.dims) && isscalar(opts.dims) && ismember(opts.dims, [2 3]))
    error('ISA:SIFTED:invalidDims', ...
        'opts.dims must be 2 or 3 (got %s).', mat2str(opts.dims));
end

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
    error('ISA:SIFTED:tooFewFeatures', ...
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
% Keep any feature with correlation at least opts.rho for any algorithm.
% Vectorised over all rows at once: applying the mask to row 1 too is
% harmless since out.selvars only ever gets set to true, never cleared,
% so it can't undo the unconditional row-1 selection above.
mask = rho >= opts.rho;
out.selvars(unique(row(mask))) = true;
out.selvars = find(out.selvars);
Xaux = X(:, out.selvars);
fprintf('[SIFTED] Keeping %d out of %d features (correlation).\n', size(Xaux,2), nfeats);
% -------------------------------------------------------------------------
nfeats = size(Xaux, 2);
if nfeats <= 1
    error('ISA:SIFTED:tooFewFeatures', ...
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
fcnwrap = @(x) costfcn(x, Xaux, Y, Ybin, out.clust, cvpart, featlabels(out.selvars), opts.dims);
% GA population fitness evaluations are parallelised at the GA level
% (UseParallel) rather than inside costfcn: a parfor over the ~10
% algorithm columns nested inside a fitness function called hundreds of
% times by the GA is nested parallelism MATLAB doesn't support (it
% silently degrades to serial) and its per-call worker-dispatch overhead
% dwarfs the tiny loop it parallelises. Parallelising GA's own
% population evaluation instead spreads real, independent work
% (one full PILOT + KNN fit per individual) across workers.
gaopts  = optimoptions('ga', ...
    'FitnessLimit',      FitnessLimit, ...
    'FunctionTolerance', FunctionTolerance, ...
    'MaxGenerations',    MaxGenerations, ...
    'MaxStallGenerations', MaxStallGenerations, ...
    'PopulationSize',    PopulationSize, ...
    'UseParallel',       nworkers > 0);
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
function y = costfcn(ind, X, Y, Ybin, clust, cvpart, featlabels, dims)
    persistent mymap
    if isempty(mymap)
        mymap = containers.Map('KeyType','char','ValueType','double');
    end
    % Internal PILOT call mirrors the canonical analytic branch, ntries=5
    % (spec §5.5), at the same dimensionality as the outer pipeline's final
    % projection. KNN uses dims+1 neighbours (D+1, generalised from 2D).
    ntries      = 5;
    analytic    = true;
    kneighbours = dims + 1;

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
        % verbose=false: this fitness function runs once per GA candidate
        % (dozens to hundreds of times per SIFTED call), so PILOT's normal
        % per-run status/summary output would flood the console.
        out = PILOT(X(:,idx), Y, featlabels(idx), ...
            struct('analytic', analytic, 'ntries', ntries, 'dims', dims, 'verbose', false));
        Z = out.Z;
        y = -Inf;
        % Plain loop, not parfor: costfcn is itself called in parallel by
        % the GA's own UseParallel population evaluation (set by the
        % caller), and nesting parfor inside that is unsupported (it
        % would silently degrade to serial) -- see the comment at the
        % gaopts definition above.
        for ii = 1:size(Y,2)
            knn = fitcknn(Z, Ybin(:,ii), 'CVPartition', cvpart, 'NumNeighbors', kneighbours);
            y = max(y, knn.kfoldLoss);
        end
        mymap(key) = y;
    end
end
% =========================================================================
function clearCache()
% Clears the persistent fitness cache in costfcn, on the client and (now
% that the GA's UseParallel dispatches costfcn calls to pool workers,
% each with its own persistent state) on every worker in the current
% parallel pool too, if one exists. Without this, a worker could return
% a stale fitness value cached under the same feature-selection bitmask
% key from a previous, unrelated SIFTED call on different data.
    clear costfcn
    if exist('gcp','file')==2
        mypool = gcp('nocreate');
        if ~isempty(mypool)
            % Clear the parent function so local-function persistent state (incl. costfcn) is reset on workers.
            wait(parfevalOnAll(mypool, @() clear('SIFTED'), 0));
        end
    end
end
% =========================================================================
