function [X, out] = SIFTED2(X, Y, Ybin, featlabels, opts)
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
% Magik numbers
Kfolds = 5; % Number of folds for the CV partition used by KNN within the GA cost function
FitnessLimit = 0; % Minimum possible for the error function
FunctionTolerance = 1e-3; % T
MaxGenerations = 100; % Maximum number of generations for the GA
MaxStallGenerations = 5; % If the GA stall for this number of generations, then exit
PopulationSize = 50; % Population size for the GA
alphaval = 0.05; % pvalue for statistical significance of the correlation
innaceptableClustering = 0.5; % A silhouette less than this value trigers a warning
acceptableClustering = 0.75; % This is a preferable silhouette value
% ---------------------------------------------------------------------

global mymap
mymap = containers.Map('KeyType','char','ValueType','double');

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
rho(isnan(rho) | (out.p>alphaval)) = 0;
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
elseif nfeats<=opts.K
    disp('-> There are less or equal features than clusters. Skipping correlation clustering selection.');
    X = Xaux;
    return;
end
% ---------------------------------------------------------------------
disp('-> Selecting features based on correlation clustering.');
state = rng;
rng('default');
out.eva = evalclusters(Xaux', 'kmeans', 'Silhouette', 'KList', 3:nfeats, ... % minimum of three features
                              'Distance', 'correlation');
disp('-> Average silhouette values for each number of clusters.')
disp([out.eva.InspectedK; out.eva.CriterionValues]);
if out.eva.CriterionValues(out.eva.InspectedK==opts.K)<innaceptableClustering
    disp(['-> The silhouette value for K=' num2str(opts.K) ...
          ' is below ' num2str(innaceptableClustering) '. You should consider increasing K.']);
    out.Ksuggested = out.eva.InspectedK(find(out.eva.CriterionValues>acceptableClustering,1));
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
disp('-> Using a GA+LookUpTable to find an optimal combination.');
% ---------------------------------------------------------------------
cvpart = cvpartition(size(Xaux,1), 'Kfold', Kfolds);
fcnwrap = @(x) costfcn(x, Xaux, Y, Ybin, out.clust, cvpart, featlabels, nworkers);
gaopts = optimoptions('ga','FitnessLimit', FitnessLimit, 'FunctionTolerance', FunctionTolerance,...
                           'MaxGenerations', MaxGenerations, 'MaxStallGenerations', MaxStallGenerations,...
                           'PopulationSize', PopulationSize); % This sets the maximum to 1000 combinations
ind = ga(fcnwrap, opts.K, [], [], [], [], ones(1,opts.K), sum(out.clust), [], 1:opts.K, gaopts);

decoder = false(1,size(Xaux,2)); % Decode the chromosome
for i=1:opts.K
    aux = find(out.clust(:,i));
    decoder(aux(ind(i))) = true;
end
out.selvars = out.selvars(decoder);
X = X(:,out.selvars);
disp(['-> Keeping ' num2str(size(X, 2)) ' out of ' num2str(nfeats) ' features (clustering).']);

end
% =========================================================================
function y = costfcn(ind,X,Y,Ybin,clust,cvpart,featlabels,nworkers)
    global mymap
    % Magik numbers
    ntries = 5; % Number of trials on PILOT
    analytic = false; % Whether the solution to PILOT is analytic or not
    kneighbours = 3; % Number of neighbours on KNN (D+1)

    idx = false(1,size(X,2));
    for i=1:length(ind)
        aux = find(clust(:,i));
        idx(aux(ind(i))) = true;
    end

    key = strrep(num2str(idx),' ','');
    if isKey(mymap,key)
        y = values(mymap,{key});
        y = y{1};
    else
        out = PILOT(X(:,idx), Y, featlabels(idx), struct('analytic', analytic, 'ntries', ntries));
        Z = out.Z;
        rng('default');
        y = -Inf;
        parfor (ii=1:size(Y,2),nworkers)
            knn = fitcknn(Z, Ybin(:,ii), 'CVPartition', cvpart, 'NumNeighbors', kneighbours);
            y = max(y, knn.kfoldLoss);
        end
        mymap(key) = y;
    end
end
% =========================================================================