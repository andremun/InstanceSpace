function [out] = matilda(X, Y, Ybin, opts)
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
%
% Input parameters:
%
% X:    A matrix of features.
% Y:    A matrix of algorithm performance.
% Ybin: A matrix with a binary measure of algorithm performance.
% opts: An structure with all the options.
%
% Output paramters:
% 
% out: An structure with all the outputs
%
% -------------------------------------------------------------------------

startProcess = tic;
% -------------------------------------------------------------------------
%
ninst = size(X,1);
nfeats = size(X,2);
nalgos = size(Y,2);
% -------------------------------------------------------------------------
%
[~,portfolio] = max(Y,[],2);
beta = sum(Ybin,2)>opts.thresholds.BETATHRESHOLD*nalgos;
% ---------------------------------------------------------------------
% Eliminate extreme outliers, i.e., any point that exceedes 5 times the
% inter quantile range, by bounding them to that value.
disp('-> Removing extreme outliers from the feature values.');
[X, out.bound] = boundOutliers(X);
% ---------------------------------------------------------------------
% Normalize the data using Box-Cox and Z-transformations
disp('-> Normalizing the data.');
[X, Y, out.norm] = autoNormalize(X, Y);
% ---------------------------------------------------------------------
% Check for diversity, i.e., we want features that have non-repeating
% values for each instance. Eliminate any that have only DIVTHRESHOLD unique values.
disp('-> Checking for feature diversity.');
[X, out.diverse] = checkDiversity(X, opts.thresholds.DIVTHRESHOLD);
disp(['-> Keeping ' num2str(size(X,2)) ' out of ' num2str(nfeats) ' features (diversity).']);
nfeats = size(X, 2);
% ---------------------------------------------------------------------
% Detect correlations between features and algorithms. Keep the top CORTHRESHOLD
% correlated features for each algorithm
disp('-> Checking for feature correlation with performance.');
[X, out.corr] = checkCorrelation(X, Y, opts.thresholds.CORTHRESHOLD);
disp(['-> Keeping ' num2str(size(X,2)) ' out of ' num2str(nfeats) ' features (correlation).']);
nfeats = size(X, 2);
% ---------------------------------------------------------------------
% Detect similar features, by clustering them according to their
% correlation. We assume that the lowest value possible is best, as this 
% will improve the projection into two dimensions. We set a hard limit of
% 10 features. The selection criteria is an average silhouete value above
% 0.65
disp('-> Selecting features based on correlation clustering.');
[X, out.clust] = clusterFeatureSelection(X, Ybin,...
                                         opts.thresholds.KDEFAULT,...
                                         opts.thresholds.SILTHRESHOLD,...
                                         opts.thresholds.NTREES);
disp(['-> Keeping ' num2str(size(X, 2)) ' out of ' num2str(nfeats) ' features (clustering).']);
nfeats = size(X, 2);
% ---------------------------------------------------------------------
% This is the final subset of features. Calculate the two dimensional
% projection using the PBLDR algorithm (Munoz et al. Mach Learn 2018)
disp('-> Finding optimum projection.');
out.pbldr = PBLDR(X, Y, opts.pbldr);
disp('-> Completed - Projection calculated. Matrix A:');
disp(out.pbldr.A);
% ---------------------------------------------------------------------
% Calculating the algorithm footprints. First step is to transform the
% data to the footprint space, and to calculate the 'space' exafootprint.
% This is also the maximum area possible for a footprint.
disp('-> Calculating the space area and density.');
spaceFootprint = findPureFootprint(out.pbldr.Z, true(ninst,1), opts.thresholds);
spaceFootprint = recalculatePerformance(spaceFootprint, out.pbldr.Z, true(ninst,1));
out.footprint.spaceArea = spaceFootprint.area;
out.footprint.spaceDensity = spaceFootprint.density;
disp(['    Area: ' num2str(out.footprint.spaceArea) ' | Density: ' num2str(out.footprint.spaceDensity)]);
% ---------------------------------------------------------------------
% This loop will calculate the footprints for good/bad instances and the
% best algorithm.
disp('-> Calculating the algorithm footprints.');
out.footprint.good = cell(1,nalgos);
out.footprint.bad = cell(1,nalgos);
out.footprint.best = cell(1,nalgos);
for i=1:nalgos
    tic;
    out.footprint.good{i} = findPureFootprint(out.pbldr.Z,...
                                              Ybin(:,i),...
                                              opts.thresholds);
    out.footprint.bad{i} = findPureFootprint( out.pbldr.Z,...
                                             ~Ybin(:,i),...
                                              opts.thresholds);
    out.footprint.best{i} = findPureFootprint(out.pbldr.Z,...
                                              portfolio==i,...
                                              opts.thresholds);
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
            out.footprint.best{i} = calculateFootprintCollisions(out.footprint.best{i}, ...
                                                                 out.footprint.best{j});
            disp(['    -> Test algorithm No. ' num2str(j) ' - Elapsed time: ' num2str(toc(startTest)) 's']);
        end
    end
    out.footprint.good{i} = calculateFootprintCollisions(out.footprint.good{i},...
                                                         out.footprint.bad{i});
    out.footprint.bad{i} = calculateFootprintCollisions(out.footprint.bad{i},...
                                                        out.footprint.good{i});
    disp(['-> Base algorithm No. ' num2str(i) ' - Elapsed time: ' num2str(toc(startBase)) 's']);
end
% -------------------------------------------------------------------------
% Calculating performance
disp('-> Calculating the footprint''s area and density. Collecting the footprint analysis results.');
out.footprint.performance = zeros(nalgos+1,10);
for i=1:nalgos
    out.footprint.best{i} = recalculatePerformance(out.footprint.best{i},...
                                                   out.pbldr.Z,...
                                                   portfolio==i);
    out.footprint.good{i} = recalculatePerformance(out.footprint.good{i},...
                                                   out.pbldr.Z,...
                                                   Ybin(:,i));
    out.footprint.bad{i} = recalculatePerformance( out.footprint.bad{i},...
                                                   out.pbldr.Z,...
                                                  ~Ybin(:,i));
    out.footprint.performance(i,:) = [calculateFootprintPerformance(out.footprint.good{i},...
                                                                    out.footprint.spaceArea,...
                                                                    out.footprint.spaceDensity), ...
                                      calculateFootprintPerformance(out.footprint.best{i},...
                                                                    out.footprint.spaceArea,...
                                                                    out.footprint.spaceDensity)];
end
% -------------------------------------------------------------------------
% Beta hard footprints. First step is to calculate them.
disp('-> Calculating beta-footprints.');
out.footprint.easy = findPureFootprint(out.pbldr.Z,  beta, opts.thresholds);
out.footprint.hard = findPureFootprint(out.pbldr.Z, ~beta, opts.thresholds);
% Remove the collisions
out.footprint.easy = calculateFootprintCollisions(out.footprint.easy,...
                                                  out.footprint.hard);
out.footprint.hard = calculateFootprintCollisions(out.footprint.hard,...
                                                  out.footprint.easy);
% Calculating performance
disp('-> Calculating the beta-footprint''s area and density. Collecting the footprint analysis results.');
out.footprint.easy = recalculatePerformance(out.footprint.easy,...
                                            out.pbldr.Z,...
                                            beta);
out.footprint.hard = recalculatePerformance(out.footprint.hard,...
                                            out.pbldr.Z,...
                                            ~beta);
out.footprint.performance(nalgos+1,:) = [calculateFootprintPerformance(out.footprint.easy,...
                                                                       out.footprint.spaceArea,...
                                                                       out.footprint.spaceDensity),...
                                         calculateFootprintPerformance(out.footprint.hard,...
                                                                       out.footprint.spaceArea,...
                                                                       out.footprint.spaceDensity)];
% ---------------------------------------------------------------------
% Making all the plots. First, plotting the features and performance as
% scatter plots.
for i=1:nfeats
    figure;
    drawScatter(X(:,i), out.pbldr.Z, ['FEATURE # ' num2str(i)]);
end

for i=1:nalgos
    figure;
    drawScatter(Y(:,i), out.pbldr.Z, ['ALGORITHM # ' num2str(i)]);
end
% Drawing the footprints for good and bad performance acording to the
% binary measure
for i=1:nalgos
    figure;
    drawGoodBadFootprint(out.footprint.good{i},...
                         out.footprint.bad{i},...
                         ['ALGORITHM # ' num2str(i)]);
end
% Drawing the footprints as portfolio.
figure;
drawPortfolioFootprint(out.footprint.best, algolabels);
% -------------------------------------------------------------------------
% 
disp(['-> Completed! Elapsed time: ' num2str(toc(startProcess)) 's']);
end
