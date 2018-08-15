function [X, out] = clusterFeatureSelection(X, Ybin, KDEFAULT, SILTHRESHOLD, NTREES)

rng('default');
nfeats = size(X,2);
nalgos = size(Ybin,2);
Kmax = min(KDEFAULT,nfeats);
out.eva = evalclusters(X', 'kmeans', 'Silhouette', 'KList', 2:Kmax, ...
                                     'Distance', 'correlation');
K = out.eva.InspectedK(find(out.eva.CriterionValues>SILTHRESHOLD,1));
out.clust = bsxfun(@eq, kmeans(X', K, 'Distance', 'correlation', ...
                                      'MaxIter', 1000, ...
                                      'Replicates', 100, ...
                                      'Options', statset('UseParallel', false), ...
                                      'OnlinePhase', 'on'), 1:K);
disp(['-> ' num2str(K) ' clusters of features detected.']);
% ---------------------------------------------------------------------
% Using these out.clusters, determine all possible combinations that take one
% feature from each out.cluster.
comb = de2bi(1:2^nfeats-1);
comb = comb(sum(comb,2)==K,:);
novalid = false(size(comb,1),1);
for i=1:K
    idx = find(out.clust(:,i));
    for j=1:length(idx)
        for k=j+1:length(idx)
            novalid(comb(:,idx(j)) & comb(:,idx(k))) = true;
        end
    end
end
comb = comb(~novalid,:)==1;
ncomb = size(comb,1);
disp(['-> ' num2str(ncomb) ' valid feature combinations.']);
% ---------------------------------------------------------------------
% Determine which combination produces the best separation while using a
% two dimensional PCA projection. The separation is defined by a Tree
% Bagger.
out.ooberr = zeros(ncomb,nalgos);
for i=1:ncomb
    tic;
    [~, score] = pca(X(:,comb(i,:)), 'NumComponents', 2); % Reduce using PCA
    for j = 1:nalgos
        rng('default');
        tree = TreeBagger(NTREES, score, Ybin(:,j), 'OOBPrediction', 'on');
        out.ooberr(i,j) = mean(Ybin(:,j)~=str2double(oobPredict(tree)));
    end
    etime = toc;
    disp(['    Combination No. ' num2str(i) ' | Elapsed Time: ' num2str(etime) ...
          ' | Average error : ' num2str(mean(out.ooberr(i,:)))]);
end
[~,best] = min(sum(out.ooberr,2));
out.selvars = comb(best,:);
X = X(:,out.selvars);
end