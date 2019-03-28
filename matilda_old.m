function [out] = matilda(featfile, algofile, opts, varargin)
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
% X:    A cell matrix of features. The first row corresponds to the feature
%       names.
% Y:    A cell matrix of algorithm performance. The first row corresponds
%       to the algorithm names.
% Ybin: A cell matrix with a binary measure of algorithm performance, where
%       '1' is good performance and '0' is bad. The first row corresponds
%       to the algorithm names.
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

% we have to change the options structure for a json file, and parse these
% into the strcuture

X = readtable(featfile);
Y = readtable(algofile);
featlabels = X.Properties.VariableNames(2:end);
algolabels = Y.Properties.VariableNames(2:end);
% -------------------------------------------------------------------------
% must compare the variable names with thoes in the JSON structure, such
% that we select the right columns to process
% -------------------------------------------------------------------------
% If we are only meant to take some observations
flag = false;
if ~isempty(varargin)
    subsetIndex = false(size(X,1),1);
    aux = readtable(varargin{1});
    subsetIndex(aux.x) = true;
    flag = true;
else
    subsetIndex = true(size(X,1),1);
end
% -------------------------------------------------------------------------
[X,idx] = sortrows(table2cell(X),1);
if flag
    subsetIndex = subsetIndex(idx);
end
Y = sortrows(table2cell(Y),1);
featfiles = X(:,1);
algofiles = Y(:,1);
X = cell2mat(X(:,2:end));
Y = cell2mat(Y(:,2:end));

isvalidY = false(length(algofiles),1);
for i=1:length(featfiles)
    isvalidY = isvalidY | strcmp(featfiles{i},algofiles);
end
algofiles = algofiles(isvalidY);
Y = Y(isvalidY,:);

isvalidX = false(length(featfiles),1);
for i=1:length(algofiles)
    isvalidX = isvalidX | strcmp(algofiles{i},featfiles);
end
X = X(isvalidX,:);
if flag
    subsetIndex = subsetIndex(isvalidX);
end
featfiles = featfiles(isvalidX);

Xbackup = X;
Ybackup = Y;
nalgos = size(Y,2);
% -------------------------------------------------------------------------
%
if opts.perf.MaxMin
    Y(isnan(Y)) = -Inf;
    [bestPerformace,portfolio] = max(Y,[],2);
    if opts.perf.AbsPerf
        Ybin = Y>=opts.perf.epsilon;
    else
        Ybin = bsxfun(@ge,Y,(1-opts.perf.epsilon).*bestPerformace); % One is good, zero is bad
    end
else
    Y(isnan(Y)) = Inf;
    [bestPerformace,portfolio] = min(Y,[],2);
    if opts.perf.AbsPerf
        Ybin = Y<=opts.perf.epsilon;
    else
        Ybin = bsxfun(@le,Y,(1+opts.perf.epsilon).*bestPerformace); 
    end
end
beta = sum(Ybin,2)>opts.general.betaThreshold*nalgos;
% ---------------------------------------------------------------------
% Eliminate extreme outliers, i.e., any point that exceedes 5 times the
% inter quantile range, by bounding them to that value.
disp('-> Removing extreme outliers from the feature values.');
[X, out.bound] = boundOutliers(X, opts.bound);
% ---------------------------------------------------------------------
% Normalize the data using Box-Cox and Z-transformations
disp('-> Auto-normalizing the data.');
[X, Y, out.norm] = autoNormalize(X, Y, opts.norm);
% ---------------------------------------------------------------------
% 
if flag
    X = X(subsetIndex,:);
    Y = Y(subsetIndex,:);
    Ybin = Ybin(subsetIndex,:);
    beta = beta(subsetIndex);
    bestPerformace = bestPerformace(subsetIndex);
    portfolio = portfolio(subsetIndex);
    Xbackup = Xbackup(subsetIndex,:);
    Ybackup = Ybackup(subsetIndex,:);
end
% -------------------------------------------------------------------------
%
ninst = size(X,1);
% ---------------------------------------------------------------------
% Check for diversity, i.e., we want features that have non-repeating
% values for each instance. Eliminate any that have only DIVTHRESHOLD unique values.
[X, out.diversity] = checkDiversity(X, opts.diversity);
featlabels = featlabels(out.diversity.selvars);
Xbackup = Xbackup(:,out.diversity.selvars);
% ---------------------------------------------------------------------
% Detect correlations between features and algorithms. Keep the top CORTHRESHOLD
% correlated features for each algorithm
[X, out.corr] = checkCorrelation(X, Y, opts.corr);
nfeats = size(X, 2);
featlabels = featlabels(out.corr.selvars);
Xbackup = Xbackup(:,out.corr.selvars);
% ---------------------------------------------------------------------
% Detect similar features, by clustering them according to their
% correlation. We assume that the lowest value possible is best, as this 
% will improve the projection into two dimensions. We set a hard limit of
% 10 features. The selection criteria is an average silhouete value above
% 0.65
disp('-> Selecting features based on correlation clustering.');
[X, out.clust] = clusterFeatureSelection(X, Ybin,...
                                         opts.clust);
disp(['-> Keeping ' num2str(size(X, 2)) ' out of ' num2str(nfeats) ' features (clustering).']);
nfeats = size(X, 2);
featlabels = featlabels(out.clust.selvars);
Xbackup = Xbackup(:,out.clust.selvars);
% ---------------------------------------------------------------------
% This is the final subset of features. Calculate the two dimensional
% projection using the PBLDR algorithm (Munoz et al. Mach Learn 2018)
disp('-> Finding optimum projection.');
out.pbldr = PBLDR(X, Y, opts.pbldr);
disp('-> Completed - Projection calculated. Matrix A:');
projectionMatrix = cell(3, nfeats+1);
projectionMatrix(1,2:end) = featlabels;
projectionMatrix(2:end,1) = {'Z_{1}','Z_{2}'};
projectionMatrix(2:end,2:end) = num2cell(out.pbldr.A);
disp(' ');
disp(projectionMatrix);
% -------------------------------------------------------------------------
% Algorithm selection. Fit a model that would separate the space into
% classes of good and bad performance. 
% out.algosel = fitoracle(out.pbldr.Z,Ybin,opts.oracle);
% ---------------------------------------------------------------------
% Calculating the algorithm footprints. First step is to transform the
% data to the footprint space, and to calculate the 'space' exafootprint.
% This is also the maximum area possible for a footprint.
disp('-> Calculating the space area and density.');
spaceFootprint = findPureFootprint(out.pbldr.Z, true(ninst,1), opts.footprint);
spaceFootprint = calculateFootprintPerformance(spaceFootprint, out.pbldr.Z, true(ninst,1));
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
                                              opts.footprint);
    out.footprint.bad{i} = findPureFootprint( out.pbldr.Z,...
                                             ~Ybin(:,i),...
                                              opts.footprint);
    out.footprint.best{i} = findPureFootprint(out.pbldr.Z,...
                                              portfolio==i,...
                                              opts.footprint);
    disp(['    -> Algorithm No. ' num2str(i) ' - Elapsed time: ' num2str(toc,'%.2f\n') 's']);
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
            disp(['    -> Test algorithm No. ' num2str(j) ' - Elapsed time: ' num2str(toc(startTest),'%.2f\n') 's']);
        end
    end
    out.footprint.good{i} = calculateFootprintCollisions(out.footprint.good{i},...
                                                         out.footprint.bad{i});
    out.footprint.bad{i} = calculateFootprintCollisions(out.footprint.bad{i},...
                                                        out.footprint.good{i});
    disp(['-> Base algorithm No. ' num2str(i) ' - Elapsed time: ' num2str(toc(startBase),'%.2f\n') 's']);
end
% -------------------------------------------------------------------------
% Calculating performance
disp('-> Calculating the footprint''s area and density.');
performanceTable = zeros(nalgos+2,10);
for i=1:nalgos
    out.footprint.best{i} = calculateFootprintPerformance(out.footprint.best{i},...
                                                   out.pbldr.Z,...
                                                   portfolio==i);
    out.footprint.good{i} = calculateFootprintPerformance(out.footprint.good{i},...
                                                   out.pbldr.Z,...
                                                   Ybin(:,i));
    out.footprint.bad{i} = calculateFootprintPerformance( out.footprint.bad{i},...
                                                   out.pbldr.Z,...
                                                  ~Ybin(:,i));
    performanceTable(i,:) = [calculateFootprintSummary(out.footprint.good{i},...
                                                       out.footprint.spaceArea,...
                                                       out.footprint.spaceDensity), ...
                             calculateFootprintSummary(out.footprint.best{i},...
                                                       out.footprint.spaceArea,...
                                                       out.footprint.spaceDensity)];
end
% -------------------------------------------------------------------------
% Beta hard footprints. First step is to calculate them.
disp('-> Calculating beta-footprints.');
out.footprint.easy = findPureFootprint(out.pbldr.Z,  beta, opts.footprint);
out.footprint.hard = findPureFootprint(out.pbldr.Z, ~beta, opts.footprint);
% Remove the collisions
out.footprint.easy = calculateFootprintCollisions(out.footprint.easy,...
                                                  out.footprint.hard);
out.footprint.hard = calculateFootprintCollisions(out.footprint.hard,...
                                                  out.footprint.easy);
% Calculating performance
disp('-> Calculating the beta-footprint''s area and density.');
out.footprint.easy = calculateFootprintPerformance(out.footprint.easy,...
                                                   out.pbldr.Z,...
                                                   beta);
out.footprint.hard = calculateFootprintPerformance(out.footprint.hard,...
                                                   out.pbldr.Z,...
                                                   ~beta);
performanceTable(end-1,6:10) = calculateFootprintSummary(out.footprint.easy,...
                                                         out.footprint.spaceArea,...
                                                         out.footprint.spaceDensity);
performanceTable(end,6:10) = calculateFootprintSummary(out.footprint.hard,...
                                                       out.footprint.spaceArea,...
                                                       out.footprint.spaceDensity);
out.footprint.performance = cell(nalgos+3,11);
out.footprint.performance(1,2:end) = {'Easy_Area',...
                                      'Easy_Area_n',...
                                      'Easy_Density',...
                                      'Easy_Density_n',...
                                      'Easy_Purity',...
                                      'Best_Area',...
                                      'Best_Area_n',...
                                      'Best_Density',...
                                      'Best_Density_n',...
                                      'Best_Purity'};
out.footprint.performance(2:end-2,1) = algolabels;
out.footprint.performance(end-1:end,1) = {'Beta-easy', 'Beta-hard'};
out.footprint.performance(2:end,2:end) = num2cell(performanceTable);
disp('-> Completed - Footprint analysis results:');
disp(' ');
disp(out.footprint.performance);
% ---------------------------------------------------------------------
% Storing the output data as a csv file
writetable(array2table(out.pbldr.Z,'VariableNames',{'z_1','z_2'},'RowNames',featfiles(subsetIndex)),...
           'coordinates.csv','WriteRowNames',true);
writetable(array2table(Xbackup,'VariableNames',featlabels,'RowNames',featfiles(subsetIndex)),...
           'feature_raw.csv','WriteRowNames',true);
writetable(array2table(X,'VariableNames',featlabels,'RowNames',featfiles(subsetIndex)),...
           'feature_process.csv','WriteRowNames',true);      
writetable(array2table(Ybackup,'VariableNames',algolabels,'RowNames',featfiles(subsetIndex)),...
           'algorithm_raw.csv','WriteRowNames',true);
writetable(array2table(Y,'VariableNames',algolabels,'RowNames',featfiles(subsetIndex)),...
           'algorithm_process.csv','WriteRowNames',true);
writetable(cell2table(out.footprint.performance(2:end,[3 5 6 8 10 11]),...
                      'VariableNames',out.footprint.performance(1,[3 5 6 8 10 11]),...
                      'RowNames',out.footprint.performance(2:end,1)),...
           'footprint_performance.csv','WriteRowNames',true);
writetable(cell2table(projectionMatrix(2:end,2:end),...
                      'VariableNames',projectionMatrix(1,2:end),...
                      'RowNames',projectionMatrix(2:end,1)),...
           'projection_matrix.csv','WriteRowNames',true);
writetable(array2table(Ybin,'VariableNames',algolabels,'RowNames',featfiles(subsetIndex)),...
           'algorithm_svm.csv','WriteRowNames',true);
% ---------------------------------------------------------------------
% Making all the plots. First, plotting the features and performance as
% scatter plots.
for i=1:nfeats
    clf;
    Xaux = (X(:,i)-min(X(:,i)))./range(X(:,i));
    aux = drawScatter(Xaux, out.pbldr.Z, strrep(strrep(featlabels{i},'_','\_'),'\_1',''));
    print(gcf,'-dpng',['scatter_' featlabels{i} '.png']);
end

Ybkup = log10(Ybackup);
Ybkup = (Ybkup-min(Ybkup(:)))./range(Ybkup(:));
for i=1:nalgos
    clf;
    Xaux = strrep(algolabels{i},'_',' ');
    Xaux = strrep(Xaux,'PERF ','');
    aux = drawScatter(Ybkup(:,i), out.pbldr.Z, Xaux);
    print(gcf,'-dpng',['scatter_' strrep(Xaux,' ','_') '.png']);
end
% Drawing the footprints for good and bad performance acording to the
% binary measure
for i=1:nalgos
    clf;
    Xaux = strrep(algolabels{i},'_',' ');
    Xaux = strrep(Xaux,'PERF ','');
    aux = drawGoodBadFootprint(out.pbldr.Z, Ybin(:,i), ...
                               out.footprint.good{i},...
                               Xaux);
    print(gcf,'-dpdf',['footprint_' strrep(Xaux,' ','_') '.pdf']);
end
% Drawing the footprints as portfolio.
clf;
aux = drawPortfolioFootprint(out.footprint.best, algolabels);
print(gcf,'-dpng','footprint_portfolio.png');

for i=1:nalgos
    clf;
    line(out.pbldr.Z(out.algosel.Yhat(:,i)==1,1), ...
         out.pbldr.Z(out.algosel.Yhat(:,i)==1,2), ...
         'LineStyle', 'none', ...
         'Marker', '.', ...
         'Color', [0.8 0.8 0.8], ...
         'MarkerFaceColor', [0.8 0.8 0.8], ...
         'MarkerSize', 6);
    line(out.pbldr.Z(out.algosel.Yhat(:,i)==2,1), ...
         out.pbldr.Z(out.algosel.Yhat(:,i)==2,2), ...
         'LineStyle', 'none', ...
         'Marker', '.', ...
         'Color', [0.0 0.0 0.0], ...
         'MarkerFaceColor', [0.0 0.0 0.0], ...
         'MarkerSize', 6);
    % handle(2) = drawFootprint(bad, [1.0 0.0 0.0]);      % Red
    xlabel('z_{1}'); ylabel('z_{2}'); title(algolabels{i});
    legend({'BAD','GOOD'}, 'Location', 'NorthEastOutside');
    set(findall(gcf,'-property','FontSize'),'FontSize',12);
    set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
    axis square;
    print(gcf,'-dpng',['footprint_' algolabels{i} '.png']);
end

% -------------------------------------------------------------------------
% 
disp(['-> Completed! Elapsed time: ' num2str(toc(startProcess)) 's']);
disp('EOF');
end
