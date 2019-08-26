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
    disp('-> Average silhouette values for each number of clusters.')
    disp([3:Kmax; out.eva.CriterionValues]);
    if all(out.eva.CriterionValues<opts.SILTHRESHOLD)
        disp('-> Silhouette threshold is too high. Using the maximum number of clusters possible.');
        K = Kmax;
    elseif all(out.eva.CriterionValues>opts.SILTHRESHOLD)
        disp('-> Silhouette threshold is too low. Consider using a higher value.');
        K = 3;
    else
        K = out.eva.InspectedK(find(out.eva.CriterionValues>opts.SILTHRESHOLD,1)); % If empty make an error that is more meaningfull
    end
    out.clust = bsxfun(@eq, kmeans(X', K, 'Distance', 'correlation', ...
                                          'MaxIter', opts.MaxIter, ...
                                          'Replicates', opts.Replicates, ...
                                          'Options', statset('UseParallel', true), ...
                                          'OnlinePhase', 'on'), 1:K);
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
    disp(['-> ' num2str(ncomb) ' valid feature combinations.']);
    % ---------------------------------------------------------------------
    % Determine which combination produces the best separation while using a
    % two dimensional PCA projection. The separation is defined by a Tree
    % Bagger.
    out.ooberr = zeros(ncomb,nalgos);
    for i=1:ncomb
        tic;
        [~, score] = pca(X(:,comb(i,:)), 'NumComponents', 2); %#ok<*IDISVAR> % Reduce using PCA
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