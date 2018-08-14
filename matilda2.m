function [out] = matilda(X,Y,Ybin)
% -------------------------------------------------------------------------
% matilda.m
% -------------------------------------------------------------------------
%
% By: Mario Andrés Muñoz Acosta
%     School of Mathematics and Statistics
%     The University of Melbourne
%     Australia
%     2018
%
% -------------------------------------------------------------------------

startProcess = tic;
% -------------------------------------------------------------------------
% Constants
DIVTHRESHOLD = 0.01; % Minimum percentage allowed of repeated values [?]
SILTHRESHOLD = 0.55; % Minimum accepted value for the average silhoute value
CORTHRESHOLD = 3;    % Top N features (by correlation) per algorithm that are selected
BETATHRESHOLD = 0.5; % Beta-easy threshold
KDEFAULT = 10;       % Default maximum number of clusters
NTREES = 50;         % Number of trees for the Random Forest (to determine highest separability in the 2-d projection)
RHO = 10;            % Density threshold
PI = 0.75;           % Purity threshold
% -------------------------------------------------------------------------
%
ninst = size(X,1);
nfeats = size(X,2);
nalgos = size(Y,2);
% -------------------------------------------------------------------------
%
[~,portfolio] = max(Y,[],2);
beta = sum(Ybin,2)>BETATHRESHOLD*nalgos;
% ---------------------------------------------------------------------
% Eliminate extreme outliers, i.e., any point that exceedes 5 times the
% inter quantile range, by bounding them to that value.
disp('-> Removing outliers from the feature values.');
X = boundOutliers(X);
% ---------------------------------------------------------------------
% Check for diversity, i.e., we want features that have non-repeating
% values for each instance. Eliminate any that have only DIVTHRESHOLD unique values.
diverse = zeros(1,nfeats);
for i=1:nfeats
    diverse(i) = length(unique(X(:,i)))./ninst;
end
hidiv = diverse>DIVTHRESHOLD;
X = X(:,hidiv);
disp(['-> Keeping ' num2str(size(X,2)) ' out of ' num2str(nfeats) ' features (diversity).']);
nfeats = size(X,2);
% ---------------------------------------------------------------------
% Normalize the data using Box-Cox and Z-transformations
disp('-> Normalizing the data.');
X = bsxfun(@minus,X,min(X,[],1))+1;
for i=1:nfeats
    X(:,i) = boxcox(X(:,i));
end
X = zscore(X);

for i=1:nalgos
    Y(:,i) = boxcox(Y(:,i));
end
Y = zscore(Y);
% ---------------------------------------------------------------------
% Detect correlations between features and algorithms. Keep the top CORTHRESHOLD
% correlated features for each algorithm
rho = corr(X,Y,'rows','pairwise');
[~,row] = sort(abs(rho),1,'descend');
hicor = unique(row(1:CORTHRESHOLD,:));
X = X(:,hicor);
disp(['-> Keeping ' num2str(size(X,2)) ' out of ' num2str(nfeats) ' features (correlation).']);
nfeats = size(X,2);
% ---------------------------------------------------------------------
% Detect similar features, by clustering them according to their
% correlation. We assume that the lowest value possible is best, as this 
% will improve the projection into two dimensions. We set a hard limit of
% 10 features. The selection criteria is an average silhouete value above
% 0.65
rng('default');
Kmax = min(KDEFAULT,nfeats);
eva = evalclusters(X', 'kmeans', 'Silhouette', 'KList', 2:Kmax, ...
                                               'Distance', 'correlation');
K = 4;%eva.InspectedK(find(eva.CriterionValues>SILTHRESHOLD,1));
clust = bsxfun(@eq, kmeans(X', K, 'Distance', 'correlation', ...
                                  'MaxIter', 1000, ...
                                  'Replicates', 100, ...
                                  'Options', statset('UseParallel', false), ...
                                  'OnlinePhase', 'on'), 1:K);
disp(['-> ' num2str(K) ' clusters of features detected.']);
% ---------------------------------------------------------------------
% Using these clusters, determine all possible combinations that take one
% feature from each cluster.
comb = de2bi(1:2^nfeats-1);
comb = comb(sum(comb,2)==K,:);
novalid = false(size(comb,1),1);
for i=1:K
    idx = find(clust(:,i));
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
ooberr = zeros(ncomb,nalgos);
for i=1:ncomb
    tic;
    [~, score] = pca(X(:,comb(i,:)), 'NumComponents', 2); % Reduce using PCA
    for j = 1:nalgos
        rng('default');
        tree = TreeBagger(NTREES, score, Ybin(:,j), 'OOBPrediction', 'on');
        ooberr(i,j) = mean(Ybin(:,j)~=str2double(oobPredict(tree)));
    end
    etime = toc;
    disp(['    Combination No. ' num2str(i) ' | Elapsed Time: ' num2str(etime) ...
          ' | Average error : ' num2str(mean(ooberr(i,:)))]);
end
[~,best] = min(sum(ooberr,2));
selvars = comb(best,:);
X = X(:,selvars);
disp(['-> Keeping ' num2str(size(X,2)) ' out of ' num2str(nfeats) ' features (clustering).']);
% ---------------------------------------------------------------------
% This is the final subset of features. Calculate the two dimensional
% projection using the PBLDR algorithm (Munoz et al. Mach Learn 2018)
disp('-> Finding optimum projection.');
opts.zflag = false;
opts.ntries = 30;
opts.fevals = 1e4;
opts.analytic = false;
[out.Z,out.A,out.B,out.C,out.error,out.R2] = PBLDR(X,Y,opts);
disp('-> Completed - Projection calculated. Matrix A:');
disp(out.A);
% ---------------------------------------------------------------------
% Calculating the algorithm footprints. First step is to transform the
% data to the footprint space, and to calculate the 'space' footprint.
% This is also the maximum area possible for a footprint.
disp('-> Calculating the space area and density.');
spaceFootprint = findPureFootprint(out.Z, false(ninst,1), RHO, PI);
spaceFootprint = recalculatePerformance(spaceFootprint, out.Z, false(ninst,1));
spaceArea = spaceFootprint.area;
spaceDensity = spaceFootprint.density;
disp(['    Area: ' num2str(spaceArea) ' | Density: ' num2str(spaceDensity)]);
% ---------------------------------------------------------------------
% Structures to hold the real and predicted footprints
disp('-> Calculating the algorithm footprints.');
solveFootprint = cell(1,nalgos);
unsolveFootprint = cell(1,nalgos);
bestFootprint = cell(1,nalgos);
% ---------------------------------------------------------------------
% This loop will calculate the footprints for solved/unsolved instances
% and the best algorithm.
for i=1:nalgos
    tic;
    solveFootprint{i} = findPureFootprint(out.Z, Ybin(:,i)==1, RHO, PI);
    unsolveFootprint{i} = findPureFootprint(out.Z, Ybin(:,i)==0, RHO, PI);
    bestFootprint{i} = findPureFootprint(out.Z, portfolio~=i, RHO, PI);
    disp(['    -> Algorithm No. ' num2str(i) ' - Elapsed time: ' num2str(toc) 's']);
end
% ---------------------------------------------------------------------
% Detecting collisions and removing them.
disp('-> Removing collisions.');
for i=1:nalgos
    startBase = tic;
    for j=1:nalgos
        if i~=j
            startTest = tic;
            bestFootprint{i} = calculateFootprintCollisions(bestFootprint{i},bestFootprint{j});
            disp(['    -> Test algorithm No. ' num2str(j) ' - Elapsed time: ' num2str(toc(startTest)) 's']);
        end
    end
    solveFootprint{i} = calculateFootprintCollisions(solveFootprint{i},unsolveFootprint{i});
    unsolveFootprint{i} = calculateFootprintCollisions(unsolveFootprint{i},solveFootprint{i});
    disp(['-> Base algorithm No. ' num2str(i) ' - Elapsed time: ' num2str(toc(startBase)) 's']);
end
% -------------------------------------------------------------------------
% Beta hard footprints. First step is to calculate them.
disp('-> Calculating beta-footprints.');
easyFootprint = findPureFootprint(out.Z,beta==1,RHO,PI);
hardFootprint = findPureFootprint(out.Z,beta==0,RHO,PI);
% The second step is to remove the collisions
easyFootprint = calculateFootprintCollisions(easyFootprint,hardFootprint);
hardFootprint = calculateFootprintCollisions(hardFootprint,easyFootprint);
% -------------------------------------------------------------------------
% Calculating performance
disp('-> Calculating the footprint''s area and density.');
for i=1:nalgos
    bestFootprint{i} = recalculatePerformance(bestFootprint{i},out.Z,portfolio~=i);
    solveFootprint{i} = recalculatePerformance(solveFootprint{i},out.Z,Ybin(:,i)==1);
    unsolveFootprint{i} = recalculatePerformance(unsolveFootprint{i},out.Z,Ybin(:,i)==0);
end
easyFootprint = recalculatePerformance(easyFootprint,out.Z,beta==1);
hardFootprint = recalculatePerformance(hardFootprint,out.Z,beta==0);
% -------------------------------------------------------------------------
% Performance of the footprints
disp('-> Collecting the performance results.');
footprintPerformance = zeros(nalgos+1,10);
% Calculating the footprint performance
for i=1:nalgos
    % Actual footprints: Solutiotion footprint
    footprintPerformance(i,:) = [calculateFootprintPerformance(solveFootprint{i},spaceArea,spaceDensity), ...
                                 calculateFootprintPerformance(bestFootprint{i},spaceArea,spaceDensity)];
end
footprintPerformance(nalgos+1,:) = [calculateFootprintPerformance(easyFootprint,spaceArea,spaceDensity),...
                                    calculateFootprintPerformance(hardFootprint,spaceArea,spaceDensity)];
% -------------------------------------------------------------------------
% 
out.footprintPerformance = footprintPerformance;
disp(['-> Completed! Elapsed time: ' num2str(toc(startProcess)) 's']);
end
% =========================================================================
% Subfunctions
% =========================================================================
function X = boundOutliers(X)
medval = median(X,1);
iqrange = iqr(X,1);
hibound = medval+5.*iqrange;
lobound = medval-5.*iqrange;
himask = bsxfun(@gt,X,hibound);
lomask = bsxfun(@lt,X,lobound);
X = X.*~(himask | lomask) + bsxfun(@times,himask,hibound) + ...
                            bsxfun(@times,lomask,lobound);
end
% -------------------------------------------------------------------------
function [Z,A,B,C,error,R2] = PBLDR(X, Y, opts)

errorfcn = @(alpha,X,n,m) mean(mean((X-(reshape(alpha((2*n)+1:end),m,2)*... % B,C
                                        reshape(alpha(1:2*n),2,n)...        % A
                                        *X(:,1:n)')').^2,1),2);

n = size(X, 2); % Number of features
Xbar = [X Y];
m = size(Xbar, 2);
if opts.zflag
    Xbar = zscore(Xbar);
    X = Xbar(:,1:n);
end
Hd = pdist(X)';

if opts.analytic
    Xbar = Xbar';
    X = X';
    [V,D] = eig(Xbar*Xbar');
    [~,idx] = sort(abs(diag(D)),'descend');
    V = V(:,idx(1:2));
    B = V(1:n,:);
    C = V(n+1:m,:)';
    Xr = X'/(X*X');
    A = V'*Xbar*Xr;
    Z = A*X;
    Xhat = [B*Z; C'*Z];
    error = sum(sum((Xbar-Xhat).^2,2));
    R2 = diag(corr(Xbar',Xhat')).^2;
else
    alpha = zeros(2*m+2*n, opts.ntries);
    eoptim = zeros(1, opts.ntries);
    perf = zeros(1, opts.ntries);

    cmaopts = bipopcmaes;
    cmaopts.StopFitness = 0;
    cmaopts.MaxRestartFunEvals = 0;
    cmaopts.MaxFunEvals  = opts.fevals;
    cmaopts.EvalParallel = 'no';
    initstr = ['2*rand(' num2str(2*m+2*n) ',1)-1'];

    for i=1:opts.ntries
        [alpha(:,i),eoptim(i)] = bipopcmaes(errorfcn, ...
                                           initstr, ...
                                           1, ...
                                           cmaopts, ...
                                           Xbar, ...
                                           n, ...
                                           m);
        A = reshape(alpha(1:2*n,i),2,n);
        Z = X*A';
        perf(i) = corr(Hd,pdist(Z)');
    end

    [~,idx] = max(perf);
    A = reshape(alpha(1:2*n,idx),2,n);
    Z = X*A';
    B = reshape(alpha((2*n)+1:end,idx),m,2);
    Xhat = Z*B';
    C = B(n+1:m,:)';
    B = B(1:n,:);
    error = sum(sum((Xbar-Xhat).^2,2));
    R2 = diag(corr(Xbar,Xhat)).^2;
end
end
% -------------------------------------------------------------------------
function footprint = findPureFootprint(x,y,rho,pi)

LOWER_PCTILE = 1;
UPPER_PCTILE = 25;

Ig        = unique(x(~y,:),'rows');   % There might be points overlapped, so eliminate them to avoid problems
numInst   = size(Ig,1);
D         = pdist(Ig)';               % Calculate the distance among all points
% Eliminate those candidates that are too small or too big
Bounds    = prctile(D,[LOWER_PCTILE UPPER_PCTILE]);
idx       = D<Bounds(1);
[row,col] = meshgrid(1:numInst);
row       = row(tril(true(numInst),-1));
col       = col(tril(true(numInst),-1));
idx       = [row(idx) col(idx)];
elim      = false(numInst,1);
for i=1:size(idx,1)
    if ~elim(idx(i,1)), elim(idx(i,2)) = true; else elim(idx(i,2)) = false; end
end
Ig          = Ig(~elim,:);
numInst     = size(Ig,1);
polygon     = delaunay(Ig);
numPolygons = size(polygon,1);
elim        = false(numPolygons,1);
for i=1:numPolygons
    elim(i) =  any(pdist(Ig(polygon(i,:),:))>Bounds(2));
end
polygon     = polygon(~elim,:);
numPolygons = size(polygon,1);
cardPoly    = zeros(numPolygons,1);
cardPolyG   = zeros(numPolygons,1);
areaPoly    = zeros(numPolygons,1);
for i=1:numPolygons
    cardPoly(i)  = sum(inpolygon(x( :,1),x( :,2),Ig(polygon(i,:),1),Ig(polygon(i,:),2)));
    cardPolyG(i) = sum(inpolygon(x(~y,1),x(~y,2),Ig(polygon(i,:),1),Ig(polygon(i,:),2)));
    areaPoly(i)  = polyarea(Ig(polygon(i,:),1),Ig(polygon(i,:),2));
end
densityPolygon = cardPoly./areaPoly;
purityPolygon  = cardPolyG./cardPoly;
idx            = (densityPolygon>=rho) & (purityPolygon>=pi);
for i=1:numPolygons
    if ~idx(i)
        sumCardPoly     = sum(cardPoly(idx));
        densityWoutPoly = sumCardPoly./sum(areaPoly(idx));
        purityWoutPoly  = sum(cardPolyG(idx))./sumCardPoly;
        if (densityWoutPoly<rho) && (purityWoutPoly<pi) %#ok<BDSCI>
            idx(i)      = true;
        end
    end
end
polygon     = polygon(idx,:);
numPolygons = size(polygon,1);
% Setting the footprint data
footprint.area    = areaPoly(idx);
footprint.density = densityPolygon(idx);
footprint.elements     = cardPoly(idx);
footprint.goodElements = cardPolyG(idx);
footprint.purity  = densityPolygon(idx);
footprint.polygon = zeros(3,2,numPolygons);
graph       = zeros(numInst);
for i=1:numPolygons
    [polyx,polyy] = poly2cw(Ig(polygon(i,:),1),Ig(polygon(i,:),2));
    footprint.polygon(:,:,i) = [polyx polyy];
    graph(polygon(i,:),polygon(i,:)) = squareform(pdist(Ig(polygon(i,:),:)));
end
graph       = sparse(graph);
distGraph   = zeros(numInst);
for i=1:numInst
    distGraph(i,:) = graphshortestpath(graph,i);
end
finish    = find(~isinf(distGraph(1,:)),1,'last');
start     = finish + 1;
footprint.pieces = 1;
while start<numInst && finish<=numInst
    finish    = start + find(~isinf(distGraph(start,start:numInst)),1,'last');
    start     = finish + 1;
    footprint.pieces = footprint.pieces + 1;
end
end
% -------------------------------------------------------------------------
function baseFootprint = calculateFootprintCollisions(baseFootprint,testFootprint)
% 
numBasePolygons = size(baseFootprint.polygon,3);
numTestPolygons = size(testFootprint.polygon,3);
totalTestArea   = zeros(numBasePolygons,numTestPolygons);
testArea        = testFootprint.area;
for i=1:numBasePolygons
    basePolygon  = baseFootprint.polygon(:,:,i);
    testPolygon1 = squeeze(testFootprint.polygon(:,1,:));
    testPolygon2 = squeeze(testFootprint.polygon(:,2,:));
    for j=1:numTestPolygons
        xx = polybool('intersection',basePolygon(:,1),basePolygon(:,2),...
                                     testPolygon1(:,j),testPolygon2(:,j));
        if ~isempty(xx), totalTestArea(i,j) = testArea(j); end
    end
end
totalTestArea = sum(totalTestArea,2);
idx           = baseFootprint.area<totalTestArea;
baseFootprint.area    = baseFootprint.area(~idx);
baseFootprint.density = baseFootprint.density(~idx);
baseFootprint.elements     = baseFootprint.elements(~idx);
baseFootprint.goodElements = baseFootprint.goodElements(~idx);
baseFootprint.purity  = baseFootprint.purity(~idx);
baseFootprint.polygon = baseFootprint.polygon(:,:,~idx);
end
% -------------------------------------------------------------------------
function footprint = recalculatePerformance(footprint,x,y)
% 
elements     = false(size(x,1),1);
goodElements = false(size(x,1),1);
for i=1:size(footprint.polygon,3)
    elements     = elements     |  inpolygon(x(:,1),x(:,2),footprint.polygon(:,1,i),footprint.polygon(:,2,i));
    goodElements = goodElements | (inpolygon(x(:,1),x(:,2),footprint.polygon(:,1,i),footprint.polygon(:,2,i)) & ~y);
end
footprint.polyArea         = footprint.area;
footprint.polyDensity      = footprint.density;
footprint.polyPurity       = footprint.purity;
footprint.polyElements     = footprint.elements;
footprint.polyGoodElements = footprint.goodElements;
footprint.elements     = elements;
footprint.goodElements = goodElements;
footprint.area         = sum(footprint.area);
footprint.density      = sum(elements)./footprint.area;
footprint.purity       = sum(goodElements)./sum(elements);
end
% -------------------------------------------------------------------------
function performance = calculateFootprintPerformance(footprint,spaceArea,spaceDensity)
% 
aux = sum(footprint.elements);
performance = zeros(1,5);
performance(1) = sum(footprint.area);
performance(2) = performance(1)/spaceArea;
performance(3) = aux./performance(1);
performance(4) = performance(3)/spaceDensity;
performance(5) = sum(footprint.goodElements)./aux;
end
% -------------------------------------------------------------------------