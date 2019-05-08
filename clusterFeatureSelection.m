function [X, out] = clusterFeatureSelection(X, Ybin, opts)

nfeats = size(X,2);
if ~opts.flag
    out.selvars = true(1,nfeats);
else

    rng('default');
    nalgos = size(Ybin,2);
    Kmax = min(opts.KDEFAULT, nfeats);
    out.eva = evalclusters(X', 'kmeans', 'Silhouette', 'KList', 3:Kmax, ... % minimum of three features
                                         'Distance', 'correlation');
    K = out.eva.InspectedK(find(out.eva.CriterionValues>opts.SILTHRESHOLD,1)); % If empty make an error that is more meaningfull
    out.clust = bsxfun(@eq, kmeans(X', K, 'Distance', 'correlation', ...
                                          'MaxIter', opts.MaxIter, ...
                                          'Replicates', opts.Replicates, ...
                                          'Options', statset('UseParallel', true), ...
                                          'OnlinePhase', 'on'), 1:K);
    disp(['-> ' num2str(K) ' clusters of features detected.']);
    % ---------------------------------------------------------------------
    % Using these out.clusters, determine all possible combinations that take one
    % feature from each out.cluster.
    strcmd = '[';
    for i=1:K
        strcmd = [strcmd 'X' num2str(i) ];
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
%     comb = de2bi(1:2^nfeats-1);
%     comb = comb(sum(comb,2)==K,:);
%     novalid = false(size(comb,1),1);
%     for i=1:K
%         idx = find(out.clust(:,i));
%         for j=1:length(idx)
%             for k=j+1:length(idx)
%                 novalid(comb(:,idx(j)) & comb(:,idx(k))) = true;
%             end
%         end
%     end
%     comb = comb(~novalid,:)==1;
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
            tree = TreeBagger(opts.NTREES, score, Ybin(:,j), 'OOBPrediction', 'on');
            out.ooberr(i,j) = mean(Ybin(:,j)~=str2double(oobPredict(tree)));
        end
        etime = toc;
        disp(['    Combination No. ' num2str(i) ' | Elapsed Time: ' num2str(etime,'%.2f\n') ...
              's | Average error : ' num2str(mean(out.ooberr(i,:)))]);
    end
    [~,best] = min(sum(out.ooberr,2));
    out.selvars = comb(best,:);
    X = X(:,out.selvars);
end
    
end