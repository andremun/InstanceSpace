function [X, out] = SIFTED(X, Y, Ybin, opts)
% -------------------------------------------------------------------------
% SIFTED.m
% -------------------------------------------------------------------------
%
% By: Mario Andres Munoz Acosta
%     School of Mathematics and Statistics
%     The University of Melbourne
%     Australia
%     2021
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

% ---------------------------------------------------------------------
nfeats = size(X,2);
if nfeats<=1
    error('-> There is only 1 feature. Stopping space construction.');
elseif nfeats<=3
    disp('-> There are 3 or less features to do selection. Skipping feature selection.')
    out.selvars = 1:nfeats;
    return;
end
% ---------------------------------------------------------------------
disp('-> Selecting features based on correlation with performance.');
[out.rho,out.p] = corr(X,Y,'rows','pairwise');
rho = out.rho;
rho(isnan(rho) | (out.p>0.05)) = 0;
[rho,row] = sort(abs(rho),1,'descend');
out.selvars = false(1,nfeats);
% Always take the most correlated feature for each algorithm
out.selvars(unique(row(1,:))) = true;
% Now take any feature that has correlation at least equal to opts.rho
for ii=2:nfeats
    out.selvars(unique(row(ii,rho(ii,:)>=opts.rho))) = true;
end
out.selvars = find(out.selvars);
Xaux = X(:,out.selvars);
disp(['-> Keeping ' num2str(size(Xaux,2)) ' out of ' num2str(nfeats) ' features (correlation).']);
% ---------------------------------------------------------------------
nfeats = size(Xaux,2);
if nfeats<=1
    error('-> There is only 1 feature. Stopping space construction.');
elseif nfeats<=3
    disp('-> There are 3 or less features to do selection. Skipping correlation clustering selection.');
    X = Xaux;
    return;
elseif nfeats<opts.K
    disp('-> There are less features than clusters. Skipping correlation clustering selection.');
    X = Xaux;
    return;
end
% ---------------------------------------------------------------------
disp('-> Selecting features based on correlation clustering.');
nalgos = size(Ybin,2);
state = rng;
rng('default');
out.eva = evalclusters(Xaux', 'kmeans', 'Silhouette', 'KList', 3:nfeats, ... % minimum of three features
                              'Distance', 'correlation');
disp('-> Average silhouette values for each number of clusters.')
disp([out.eva.InspectedK; out.eva.CriterionValues]);
if out.eva.CriterionValues(out.eva.InspectedK==opts.K)<0.5
    disp(['-> The silhouette value for K=' num2str(opts.K) ...
          ' is below 0.5. You should consider increasing K.']);
    out.Ksuggested = out.eva.InspectedK(find(out.eva.CriterionValues>0.75,1));
    if ~isempty(out.Ksuggested)
        disp(['-> A suggested value of K is ' num2str(out.Ksuggested)]);
    end
end
% ---------------------------------------------------------------------
rng('default');
out.clust = bsxfun(@eq, kmeans(Xaux', opts.K, 'Distance', 'correlation', ...
                                              'MaxIter', opts.MaxIter, ...
                                              'Replicates', opts.Replicates, ...
                                              'Options', statset('UseParallel', nworkers~=0), ...
                                              'OnlinePhase', 'on'), 1:opts.K);
rng(state);
disp(['-> Constructing ' num2str(opts.K) ' clusters of features.']);
% ---------------------------------------------------------------------
% Using these out.clusters, determine all possible combinations that take one
% feature from each out.cluster.
strcmd = '[';
for i=1:opts.K
    strcmd = [strcmd 'X' num2str(i) ]; %#ok<*AGROW>
    if i<opts.K
        strcmd = [strcmd ','];
    else
        strcmd = [strcmd '] = ndgrid('];
    end
end
for i=1:opts.K
    strcmd = [strcmd 'find(out.clust(:,' num2str(i) '))'];
    if i<opts.K
        strcmd = [strcmd ','];
    else
        strcmd = [strcmd ');'];
    end
end
eval(strcmd);
strcmd = 'comb = [';
for i=1:opts.K
    strcmd = [strcmd 'X' num2str(i) '(:)'];
    if i<opts.K
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
    fcnwrap = @(idx) fcnforga(idx,Xaux,Ybin,opts.NTREES,out.clust,nworkers);
    gaopts = optimoptions('ga','FitnessLimit',0,'FunctionTolerance',1e-3,...
                          'MaxGenerations',100,'MaxStallGenerations',5,...
                          'PopulationSize',50); % This sets the maximum to 1000 combinations
    ind = ga(fcnwrap,opts.K,[],[],[],[],ones(1,opts.K),sum(out.clust),[],1:opts.K,gaopts);
    decoder = false(1,size(Xaux,2)); % Decode the chromosome
    for i=1:opts.K
        aux = find(out.clust(:,i));
        decoder(aux(ind(i))) = true;
    end
    out.selvars = out.selvars(decoder);
elseif ncomb==1
    disp('-> There is one valid combination. It will be considered the optimal one.');
    % out.selvars = 1:nfeats;
else
    disp('-> There are less than 1000 valid combinations. Using brute-force to find an optimal one.');
    out.ooberr = zeros(ncomb,nalgos);
    for i=1:ncomb
        tic;
        out.ooberr(i,:) = costfcn(comb(i,:),Xaux,Ybin,opts.NTREES,nworkers);
        etime = toc;
        disp(['    -> Combination No. ' num2str(i) ' | Elapsed Time: ' num2str(etime,'%.2f\n') ...
              's | Average error : ' num2str(mean(out.ooberr(i,:)))]);
        tic;
    end
    [~,best] = min(sum(out.ooberr,2));
    out.selvars = sort(out.selvars(comb(best,:)));
end
X = X(:,out.selvars);
disp(['-> Keeping ' num2str(size(X, 2)) ' out of ' num2str(nfeats) ' features (clustering).']);

end
% =========================================================================
function ooberr = costfcn(comb,X,Ybin,ntrees,nworkers)

[~, score] = pca(X(:,comb), 'NumComponents', 2); %#ok<*IDISVAR> % Reduce using PCA
nalgos = size(Ybin,2);
ooberr = zeros(1,nalgos);
for j = 1:nalgos
    state = rng;
    rng('default');
    tree = TreeBagger(ntrees, score, Ybin(:,j), 'OOBPrediction', 'on',...
                                                'Options', statset('UseParallel', nworkers~=0));
    ooberr(j) = mean(Ybin(:,j)~=str2double(oobPredict(tree)));
    rng(state);
end

end
% =========================================================================
function sumerr = fcnforga(idx,X,Ybin,ntrees,clust,nworkers)
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
    ooberr = costfcn(ccomb,X,Ybin,ntrees,nworkers);
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