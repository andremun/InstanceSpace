function [X, out] = fitInstanceSpace(X, Y, Ybin, opts)

nfeats = size(X,2);
rng('default');
Kmax = min(opts.clust.KDEFAULT, nfeats);
out.clust.eva = evalclusters(X', 'kmeans', 'Silhouette', 'KList', 2:Kmax, ...
                                 'Distance', 'correlation');
K = out.clust.eva.InspectedK(find(out.clust.eva.CriterionValues>opts.clust.SILTHRESHOLD,1)); % If empty make an error that is more meaningfull
out.clust.cidx = bsxfun(@eq, kmeans(X', K, 'Distance', 'correlation', ...
                                           'MaxIter', opts.clust.MaxIter, ...
                                           'Replicates', opts.clust.Replicates, ...
                                           'Options', statset('UseParallel', true), ...
                                           'OnlinePhase', 'on'), 1:K);
disp(['-> ' num2str(K) ' clusters of features detected.']);
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
    strcmd = [strcmd 'find(out.clust.cidx(:,' num2str(i) '))'];
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
% PBLDR projection. The separation is defined by a crossvalidated SVM.
ccvmcr = Inf;
disp('-> Testing each combination');
for i=1:ncomb
    tic;
    disp('-> Finding optimum projection.');
    aux.pbldr = PBLDR(X(:,comb(i,:)), Y, opts.pbldr);
    disp('-> Testing space separation.');
    aux.oracle = fitoracle(aux.pbldr.Z, Ybin, opts.oracle);
    aux.cvmcr = mean(min(aux.oracle.cvmcr,[],1));
    if aux.cvmcr<ccvmcr
        disp('-> Best combination so far found.');
        ccvmcr = aux.cvmcr;
        out.selvars = comb(i,:);
        out.pbldr = aux.pbldr;
        out.oracle = aux.oracle;
    end
    etime = toc;
    disp(['    Combination No. ' num2str(i) ' | Elapsed Time: ' num2str(etime,'%.2f\n') ...
          's | Average error : ' num2str(aux.cvmcr)]);
end
X = X(:,out.selvars);
end