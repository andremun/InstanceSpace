function [X, out] = clusterFeatureSelection(X, Ybin, opts)
% -------------------------------------------------------------------------
% clusterFeatureSelection.m
% -------------------------------------------------------------------------
%
% By: Mario Andres Munoz Acosta
%     School of Mathematics and Statistics
%     The University of Melbourne
%     Australia
%     2019
%
% -------------------------------------------------------------------------

global comb gacostvals

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

nfeats = size(X,2);
disp('-> Selecting features based on correlation clustering.');
nalgos = size(Ybin,2);
Kmax = min(opts.KDEFAULT, nfeats);

state = rng;
rng('default');
out.eva = evalclusters(X', 'kmeans', 'Silhouette', 'KList', 3:Kmax, ... % minimum of three features
                           'Distance', 'correlation');
rng(state);
disp('-> Average silhouette values for each number of clusters.')
disp([3:Kmax; out.eva.CriterionValues]);

if all(out.eva.CriterionValues<opts.SILTHRESHOLD)
    disp('-> Silhouette threshold is too high. Consider using a lower value.');
    disp('-> Using the maximum number of clusters possible.');
    K = Kmax;
elseif all(out.eva.CriterionValues>opts.SILTHRESHOLD)
    disp('-> Silhouette threshold is too low. Consider using a higher value.');
    disp('-> Using the minimum number of clusters possible.');
    K = 3;
else
    K = out.eva.InspectedK(find(out.eva.CriterionValues>opts.SILTHRESHOLD,1)); % If empty make an error that is more meaningfull
end

state = rng;
rng('default');
out.clust = bsxfun(@eq, kmeans(X', K, 'Distance', 'correlation', ...
                                      'MaxIter', opts.MaxIter, ...
                                      'Replicates', opts.Replicates, ...
                                      'Options', statset('UseParallel', nworkers~=0), ...
                                      'OnlinePhase', 'on'), 1:K);
rng(state);
disp(['-> Constructing ' num2str(K) ' clusters of features.']);
% ---------------------------------------------------------------------
% Using these out.clusters, determine all possible combinations that take one
% feature from each out.cluster.
strcmd = '[';
for i=1:K
    strcmd = [strcmd 'X' num2str(i) ]; %#ok<*AGROW>
    if i<K
        strcmd = [strcmd ','];
    else
        strcmd = [strcmd '] = ndgrid('];
    end
end
for i=1:K
    strcmd = [strcmd 'find(out.clust(:,' num2str(i) '))'];
    if i<K
        strcmd = [strcmd ','];
    else
        strcmd = [strcmd ');'];
    end
end
eval(strcmd);
strcmd = 'comb = [';
for i=1:K
    strcmd = [strcmd 'X' num2str(i) '(:)'];
    if i<K
        strcmd = [strcmd ','];
    else
        strcmd = [strcmd '];'];
    end
end
eval(strcmd);

ncomb = size(comb,1); %#ok<*NODEF>
comb = sort(comb,2);
disp(['-> ' num2str(ncomb) ' valid feature combinations.']);

maxcomb = 1000;
% ---------------------------------------------------------------------
% Determine which combination produces the best separation while using a
% two dimensional PCA projection. The separation is defined by a Tree
% Bagger.
if ncomb>maxcomb
    disp('-> There are over 1000 valid combinations. Using a GA+LookUpTable to find an optimal one.');
    gacostvals = NaN.*ones(ncomb,1);
    % gacostvals = spalloc(1,bi2de(ones(1,nfeats)),ncomb);
    fcnwrap = @(idx) fcnforga(idx,X,Ybin,opts.NTREES,out.clust);
    gaopts = optimoptions('ga','FitnessLimit',0,'FunctionTolerance',1e-3,...
                          'MaxGenerations',100,'MaxStallGenerations',5,...
                          'PopulationSize',50); % This sets the maximum to 1000 combinations
    ind = ga(fcnwrap,K,[],[],[],[],ones(1,K),sum(out.clust),[],1:K,gaopts);
    out.selvars = false(1,size(X,2)); % Decode the chromosome
    for i=1:K
        aux = find(out.clust(:,i));
        out.selvars(aux(ind(i))) = true;
    end
    out.selvars = find(out.selvars);
elseif ncomb==1
    disp('-> There is one valid combination. It will be considered the optimal one.');
    out.selvars = 1:nfeats;
else
    disp('-> There are less than 1000 valid combinations. Using brute-force to find an optimal one.');
    out.ooberr = zeros(ncomb,nalgos);
    for i=1:ncomb
        tic;
        out.ooberr(i,:) = costfcn(comb(i,:),X,Ybin,opts.NTREES);
        etime = toc;
        disp(['    -> Combination No. ' num2str(i) ' | Elapsed Time: ' num2str(etime,'%.2f\n') ...
              's | Average error : ' num2str(mean(out.ooberr(i,:)))]);
        tic;
    end
    [~,best] = min(sum(out.ooberr,2));
    out.selvars = sort(comb(best,:));
end
X = X(:,out.selvars);
disp(['-> Keeping ' num2str(size(X, 2)) ' out of ' num2str(nfeats) ' features (clustering).']);

end
% =========================================================================
function ooberr = costfcn(comb,X,Ybin,ntrees)

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

[~, score] = pca(X(:,comb), 'NumComponents', 2); %#ok<*IDISVAR> % Reduce using PCA
nalgos = size(Ybin,2);
ooberr = zeros(1,nalgos);
for j = 1:nalgos
    state = rng;
    rng('default');
    tree = TreeBagger(ntrees, score, Ybin(:,j), 'OOBPrediction', 'on', ...
                      'Options', statset('UseParallel', nworkers~=0));
    ooberr(j) = mean(Ybin(:,j)~=str2double(oobPredict(tree)));
    rng(state);
end

end
% =========================================================================
function sumerr = fcnforga(idx,X,Ybin,ntrees,clust)
global gacostvals comb

tic;
% Decode the chromosome into a binary string representing the selected
% features
ccomb = false(1,size(X,2));
for i=1:length(idx)
    aux = find(clust(:,i));
    ccomb(aux(idx(i))) = true;
end
% Calculate the cost function. Use the lookup table to reduce the amount of
% computation.
ind = find(all(comb==find(ccomb),2));
if isnan(gacostvals(ind))
    ooberr = costfcn(ccomb,X,Ybin,ntrees);
    sumerr = sum(ooberr);
    gacostvals(ind) = sumerr;
else
    sumerr = gacostvals(ind);
end

etime = toc;
disp(['    -> Combination No. ' num2str(ind) ' | Elapsed Time: ' num2str(etime,'%.2f\n') ...
      's | Average error : ' num2str(sumerr./size(Ybin,2))]);

end
% =========================================================================