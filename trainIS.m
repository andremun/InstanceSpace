function [model] = trainIS(rootdir)
% -------------------------------------------------------------------------
% trainIS.m
% -------------------------------------------------------------------------
%
% By: Mario Andres Munoz Acosta
%     School of Mathematics and Statistics
%     The University of Melbourne
%     Australia
%     2019
%
% -------------------------------------------------------------------------

startProcess = tic;
printdisclaim('trainIS.m');
% -------------------------------------------------------------------------
% Collect all the data from the files
datafile = [rootdir 'metadata.csv'];
optsfile = [rootdir 'options.json'];
if ~isfile(datafile) || ~isfile(optsfile)
    error('Please place the datafiles in the directory specified.');
end
opts = jsondecode(fileread(optsfile));
disp('-------------------------------------------------------------------------');
disp('-> Listing options to be used:');
optfields = fieldnames(opts);
for i = 1:length(optfields)
    disp(optfields{i});
    disp(opts.(optfields{i}));
end
disp('-------------------------------------------------------------------------');
disp('-> Loading the data.');
Xbar = readtable(datafile);
varlabels = Xbar.Properties.VariableNames;
isname = strcmpi(varlabels,'instances');
isfeat = strncmpi(varlabels,'feature_',8);
isalgo = strncmpi(varlabels,'algo_',5);
issource = strcmpi(varlabels,'source');
instlabels = Xbar{:,isname};
if isnumeric(instlabels)
    instlabels = num2cell(instlabels);
    instlabels = cellfun(@(x) num2str(x),instlabels,'UniformOutput',false);
end
if any(issource)
    S = categorical(Xbar{:,issource});
end
X = Xbar{:,isfeat};
Y = Xbar{:,isalgo};
% -------------------------------------------------------------------------
% Giving the oportunity to pick and choose which features/algorithms to
% work with
featlabels = varlabels(isfeat);
if isfield(opts,'selvars') && isfield(opts.selvars,'feats')
    disp('-------------------------------------------------------------------------');
    msg = '-> Using the following features: ';
    isselfeat = false(1,length(featlabels));
    for i=1:length(opts.selvars.feats)
        isselfeat = isselfeat | strcmp(featlabels,opts.selvars.feats{i});
        msg = [msg opts.selvars.feats{i} ' ']; %#ok<AGROW>
    end
    disp(msg);
    X = X(:,isselfeat);
    featlabels = featlabels(isselfeat);
end

algolabels = varlabels(isalgo);
if isfield(opts,'selvars') && isfield(opts.selvars,'algos')
    disp('-------------------------------------------------------------------------');
    msg = '-> Using the following algorithms: ';
    isselalgo = false(1,length(algolabels));
    for i=1:length(opts.selvars.algos)
        isselalgo = isselalgo | strcmp(algolabels,opts.selvars.algos{i});
        msg = [msg opts.selvars.algos{i} ' ']; %#ok<AGROW>
    end
    disp(msg);
    Y = Y(:,isselalgo);
    algolabels = algolabels(isselalgo);
end
nalgos = size(Y,2);
% -------------------------------------------------------------------------
% Storing the raw data for further processing, e.g., graphs
Xraw = X;
Yraw = Y;
% -------------------------------------------------------------------------
% Removing the template data such that it can be used in the labels of
% graphs and figures.
featlabels = strrep(featlabels,'feature_','');
algolabels = strrep(algolabels,'algo_','');
% -------------------------------------------------------------------------
% Determine whether the performance of an algorithm is a cost measure to
% be minimized or a profit measure to be maximized. Moreover, determine
% whether we are using an absolute threshold as good peformance (the
% algorithm has a performance better than the threshold) or a relative
% performance (the algorithm has a performance that is similar that the
% best algorithm minus a percentage).
disp('-------------------------------------------------------------------------');
disp('-> Calculating the binary measure of performance');
msg = '-> An algorithm is good if its performace is ';
if opts.perf.MaxMin
    Y(isnan(Y)) = -Inf;
    [bestPerformace,portfolio] = max(Y,[],2);
    if opts.perf.AbsPerf
        Ybin = Y>=opts.perf.epsilon;
        msg = [msg 'higher than ' num2str(opts.perf.epsilon)];
    else
        Ybin = bsxfun(@ge,Y,(1-opts.perf.epsilon).*bestPerformace); % One is good, zero is bad
        msg = [msg 'within ' num2str(round(100.*opts.perf.epsilon)) '% of the best.'];
    end
else
    Y(isnan(Y)) = Inf;
    [bestPerformace,portfolio] = min(Y,[],2);
    if opts.perf.AbsPerf
        Ybin = Y<=opts.perf.epsilon;
        msg = [msg 'less than ' num2str(opts.perf.epsilon)];
    else
        Ybin = bsxfun(@le,Y,(1+opts.perf.epsilon).*bestPerformace);
        msg = [msg 'within ' num2str(round(100.*opts.perf.epsilon)) '% of the best.'];
    end
end
portfolio(sum(bsxfun(@eq,Y,bestPerformace),2)>1) = 0;
disp(msg);
beta = sum(Ybin,2)>opts.general.betaThreshold*nalgos;
% ---------------------------------------------------------------------
% Automated pre-processing
if opts.auto.preproc
    disp('-------------------------------------------------------------------------');
    disp('-> Auto-pre-processing.');
    % Eliminate extreme outliers, i.e., any point that exceedes 5 times the
    % inter quantile range, by bounding them to that value.
    [X, model.bound] = boundOutliers(X, opts.bound);
    % Normalize the data using Box-Cox and Z-transformations
    [X, Y, model.norm] = autoNormalize(X, Y, opts.norm);
end
% -------------------------------------------------------------------------
% If we are only meant to take some observations
disp('-------------------------------------------------------------------------');
ninst = size(X,1);
fractional = opts.selvars.smallscaleflag && isfloat(opts.selvars.smallscale);
fileindexed = opts.selvars.fileidxflag && isfield(opts,'selvars') && isfield(opts.selvars,'instances') && isfile(opts.selvars.fileidx);
if fractional
    disp(['-> Creating a small scale experiment for validation. Percentage of subset: ' ...
        num2str(round(100.*opts.selvars.smallscale,2)) '%']);
    aux = cvpartition(ninst,'HoldOut',opts.selvars.smallscale);
    subsetIndex = aux.test;
elseif fileindexed
    disp('-> Using a subset of the instances.');
    subsetIndex = false(size(X,1),1);
    aux = table2array(readtable(opts.selvars.instances));
    aux(aux>ninst) = [];
    subsetIndex(aux) = true;
else
    disp('-> Using the complete set of the instances.');
    subsetIndex = true(ninst,1);
end

if fileindexed || fractional
    X = X(subsetIndex,:);
    Y = Y(subsetIndex,:);
    Ybin = Ybin(subsetIndex,:);
    beta = beta(subsetIndex);
    bestPerformace = bestPerformace(subsetIndex); 
    portfolio = portfolio(subsetIndex);
end
[ninst,nfeats] = size(X);
% ---------------------------------------------------------------------
% Automated feature selection.
% Keep track of the features that have been removed so we can use them
% later
model.featsel.idx = 1:nfeats;
if opts.auto.featsel
    disp('-------------------------------------------------------------------------');
    disp('-> Auto-feature selection.');
    % Check for diversity, i.e., we want features that have non-repeating
    % values for each instance. Eliminate any that have only DIVTHRESHOLD
    % unique values.
    [X, model.diversity] = checkDiversity(X, opts.diversity);
    featlabels = featlabels(model.diversity.selvars);
    model.featsel.idx = model.featsel.idx(model.diversity.selvars);
    % Detect correlations between features and algorithms. Keep the top
    % CORTHRESHOLD correlated features for each algorithm
    [X, model.corr] = checkCorrelation(X, Y, opts.corr);
    featlabels = featlabels(model.corr.selvars);
    model.featsel.idx = model.featsel.idx(model.corr.selvars);
    % Detect similar features, by clustering them according to their
    % correlation. We assume that the lowest value possible is best, as
    % this will improve the projection into two dimensions. We set a hard
    % limit of 10 features. The selection criteria is an average silhouete
    % value above 0.65
    [X, model.clust] = clusterFeatureSelection(X, Ybin, opts.clust);
    nfeats = size(X, 2);
    featlabels = featlabels(model.clust.selvars);
    model.featsel.idx = model.featsel.idx(model.clust.selvars);
end
% ---------------------------------------------------------------------
% This is the final subset of features. Calculate the two dimensional
% projection using the PBLDR algorithm (Munoz et al. Mach Learn 2018)
disp('-------------------------------------------------------------------------');
disp('-> Finding optimum projection.');
model.pbldr = PBLDR(X, Y, opts.pbldr);
disp('-> Completed - Projection calculated. Matrix A:');
projectionMatrix = cell(3, nfeats+1);
projectionMatrix(1,2:end) = featlabels;
projectionMatrix(2:end,1) = {'Z_{1}','Z_{2}'};
projectionMatrix(2:end,2:end) = num2cell(round(model.pbldr.A,4));
disp(' ');
disp(projectionMatrix);
% ---------------------------------------------------------------------
% Finding the empirical bounds based on the ranges of the features and the
% correlations of the different edges.
disp('-------------------------------------------------------------------------');
disp('-> Finding empirical bounds.');
model.sbound = findSpaceBounds(X,model.pbldr.A,opts.sbound);
% % -------------------------------------------------------------------------
% % Algorithm selection. Fit a model that would separate the space into
% % classes of good and bad performance. 
% disp('-------------------------------------------------------------------------');
% disp('-> Fitting SVM prediction models. This may take a while...');
% model.algosel = fitoracle(model.pbldr.Z, Ybin, opts.oracle);
% % svmselections = bsxfun(@eq,model.algosel.psel,unique(model.algosel.psel)');
% svmselections = bsxfun(@eq,model.algosel.psel,1:nalgos);
% Yaux = Yraw(subsetIndex,:);
% Yaux(~svmselections) = NaN;
% svmTable = cell(8,nalgos+3);
% svmTable{1,1} = ' ';
% svmTable(1,2:end-2) = algolabels;
% svmTable(1,end-1:end) = {'Oracle','Selector'};
% svmTable(2:8,1) = {'Avg. Perf. all instances';
%                    'Std. Perf. all instances';
%                    'Avg. Perf. selected instances';
%                    'Std. Perf. selected instances';
%                    'CV model error';
%                    'C';
%                    'Gamma'};
% svmTable(2,2:end) = num2cell([mean(Yraw(subsetIndex,:)) mean(bestPerformace) nanmean(Yaux(:))]);
% svmTable(3,2:end) = num2cell([std(Yraw(subsetIndex,:)) std(bestPerformace) nanstd(Yaux(:))]);
% svmTable(4,2:end-2) = num2cell(nanmean(Yaux));
% svmTable(5,2:end-2) = num2cell(nanstd(Yaux));
% svmTable(6,2:end-2) = num2cell(round(100.*model.algosel.modelerr,1));
% svmTable(7,2:end-2) = num2cell(model.algosel.svmparams(:,1));
% svmTable(8,2:end-2) = num2cell(model.algosel.svmparams(:,2));
% disp('-> Completed! Performance of the models:');
% disp(' ');
% disp(svmTable);
% ---------------------------------------------------------------------
% Calculating the algorithm footprints. First step is to transform the
% data to the footprint space, and to calculate the 'space' exafootprint.
% This is also the maximum area possible for a footprint.
disp('-------------------------------------------------------------------------');
disp('-> Calculating the space area and density.');
spaceFootprint = findPureFootprint(model.pbldr.Z, true(ninst,1), opts.footprint);
spaceFootprint = calculateFootprintPerformance(spaceFootprint, model.pbldr.Z, true(ninst,1));
model.footprint.spaceArea = spaceFootprint.area;
model.footprint.spaceDensity = spaceFootprint.density;
disp(['    Space area: ' num2str(model.footprint.spaceArea) ' | Space density: ' num2str(model.footprint.spaceDensity)]);
% ---------------------------------------------------------------------
% This loop will calculate the footprints for good/bad instances and the
% best algorithm.
disp('-------------------------------------------------------------------------');
disp('-> Calculating the algorithm footprints.');
model.footprint.good = cell(1,nalgos);
model.footprint.bad = cell(1,nalgos);
model.footprint.best = cell(1,nalgos);
% if opts.footprint.usesim
%     % Use the SVM prediction data to calculate the footprints
%     disp('-> Using prediction data.');
%     Yfoot = model.algosel.Yhat;
%     Pfoot = model.algosel.psel;
%     Bfoot = sum(model.algosel.Yhat,2)>opts.general.betaThreshold*nalgos;
% else
    % Use the actual data to calculate the footprints
    disp('-> Using experimental data.');
    Yfoot = Ybin;
    Pfoot = portfolio;
    Bfoot = beta;
% end

for i=1:nalgos
    tic;
    disp(['    -> Good performance footprint for ''' algolabels{i} '''']);
    model.footprint.good{i} = findPureFootprint(model.pbldr.Z,...
                                                Yfoot(:,i),...
                                                opts.footprint);
    disp(['    -> Bad performance footprint for ''' algolabels{i} '''']);
    model.footprint.bad{i} = findPureFootprint( model.pbldr.Z,...
                                               ~Yfoot(:,i),...
                                               opts.footprint);
    disp(['    -> Best performance footprint for ''' algolabels{i} '''']);
    model.footprint.best{i} = findPureFootprint(model.pbldr.Z,...
                                                Pfoot==i,...
                                                opts.footprint);
    disp(['    -> Algorithm ''' algolabels{i} ''' completed - Elapsed time: ' num2str(toc,'%.2f\n') 's']);
end
% ---------------------------------------------------------------------
% Detecting collisions and removing them.
disp('-------------------------------------------------------------------------');
disp('-> Detecting and removing contradicting sections of the footprints.');
for i=1:nalgos
    disp(['  -> Base algorithm ''' algolabels{i} '''']);
    startBase = tic;
    for j=i+1:nalgos
        disp(['    -> Comparing ''' algolabels{i} ''' with ''' algolabels{j} '''']);
        startTest = tic;
        [model.footprint.best{i},...
         model.footprint.best{j}] = calculateFootprintCollisionsDual(model.footprint.best{i}, ...
                                                                     model.footprint.best{j});
        
        disp(['    -> Test algorithm ''' algolabels{j} ...
              ''' completed - Elapsed time: ' num2str(toc(startTest),'%.2f\n') 's']);
    end
    disp(['   -> Comparing good and bad performance areas for ''' algolabels{i} '''']);
    [model.footprint.good{i},...
     model.footprint.bad{i}] = calculateFootprintCollisionsDual(model.footprint.good{i},...
                                                                model.footprint.bad{i});
    disp(['  -> Base algorithm ''' algolabels{i} ...
          ''' completed - Elapsed time: ' num2str(toc(startBase),'%.2f\n') 's']);
end
% -------------------------------------------------------------------------
% Calculating performance
disp('-------------------------------------------------------------------------');
disp('-> Calculating the footprint''s area and density.');
performanceTable = zeros(nalgos+2,10);
for i=1:nalgos
    model.footprint.best{i} = calculateFootprintPerformance(model.footprint.best{i},...
                                                            model.pbldr.Z,...
                                                            Pfoot==i);
    model.footprint.good{i} = calculateFootprintPerformance(model.footprint.good{i},...
                                                            model.pbldr.Z,...
                                                            Yfoot(:,i));
    model.footprint.bad{i} = calculateFootprintPerformance( model.footprint.bad{i},...
                                                            model.pbldr.Z,...
                                                           ~Yfoot(:,i));
    performanceTable(i,:) = [calculateFootprintSummary(model.footprint.good{i},...
                                                       model.footprint.spaceArea,...
                                                       model.footprint.spaceDensity), ...
                             calculateFootprintSummary(model.footprint.best{i},...
                                                       model.footprint.spaceArea,...
                                                       model.footprint.spaceDensity)];
    disp(['    -> Algorithm ' algolabels{i} ' completed - Elapsed time: ' num2str(toc,'%.2f\n') 's']);
end
% -------------------------------------------------------------------------
% Beta hard footprints. First step is to calculate them.
disp('-------------------------------------------------------------------------');
disp('-> Calculating beta-footprints.');
model.footprint.easy = findPureFootprint(model.pbldr.Z,  Bfoot, opts.footprint);
model.footprint.hard = findPureFootprint(model.pbldr.Z, ~Bfoot, opts.footprint);
% Remove the collisions
[model.footprint.easy,...
 model.footprint.hard] = calculateFootprintCollisionsDual(model.footprint.easy,...
                                                          model.footprint.hard);
% Calculating performance
disp('-> Calculating the beta-footprint''s area and density.');
model.footprint.easy = calculateFootprintPerformance(model.footprint.easy,...
                                                     model.pbldr.Z,...
                                                     Bfoot);
model.footprint.hard = calculateFootprintPerformance( model.footprint.hard,...
                                                      model.pbldr.Z,...
                                                     ~Bfoot);
performanceTable(end-1,6:10) = calculateFootprintSummary(model.footprint.easy,...
                                                         model.footprint.spaceArea,...
                                                         model.footprint.spaceDensity);
performanceTable(end,6:10) = calculateFootprintSummary(model.footprint.hard,...
                                                       model.footprint.spaceArea,...
                                                       model.footprint.spaceDensity);
model.footprint.performance = cell(nalgos+3,11);
disp('-------------------------------------------------------------------------');
model.footprint.performance(1,2:end) = {'Area_Good',...
                                        'Area_Good_Normalized',...
                                        'Density_Good',...
                                        'Density_Good_Normalized',...
                                        'Purity_Good',...
                                        'Area_Best',...
                                        'Area_Best_Normalized',...
                                        'Density_Best',...
                                        'Density_Best_Normalized',...
                                        'Purity_Best'};
model.footprint.performance(2:end-2,1) = algolabels;
model.footprint.performance(end-1:end,1) = {'beta-easy', 'beta-hard'};
model.footprint.performance(2:end,2:end) = num2cell(round(performanceTable,3));
disp('-> Completed - Footprint analysis results:');
disp(' ');
disp(model.footprint.performance);
% -------------------------------------------------------------------------
% Algorithm selection. Fit a model that would separate the space into
% classes of good and bad performance. 
disp('-------------------------------------------------------------------------');
disp('-> Fitting SVM prediction models. This may take a while...');
model.algosel = fitoracle(model.pbldr.Z, Ybin, ...
                          cell2mat(model.footprint.performance(2:nalgos+1,3)), ...
                          opts.oracle);
% svmselections = bsxfun(@eq,model.algosel.psel,unique(model.algosel.psel)');
svmselections = bsxfun(@eq,model.algosel.psel,1:nalgos);
Yaux = Yraw(subsetIndex,:);
Yaux(~svmselections) = NaN;
svmTable = cell(nalgos+3, 8);
svmTable{1,1} = ' ';
svmTable(2:end-2, 1)= algolabels;
svmTable(end-1:end, 1)={'Oracle','Selector'};
svmTable(1, 2:8)= {'Average_Performance_all_instances';
                   'Standard_Deviation_Performance_all_instances';
                   'Average_Performance_selected_instances';
                   'Standard_Deviation_Performance_selected_instances';
                   'CV_model_accuracy';
                   'C';
                   'Gamma'};
svmTable(2:end,2) = num2cell([mean(Yraw(subsetIndex,:)) mean(bestPerformace) nanmean(Yaux(:))]);
svmTable(2:end,3) = num2cell([std(Yraw(subsetIndex,:)) std(bestPerformace) nanstd(Yaux(:))]);
svmTable(2:end-2,4) = num2cell(nanmean(Yaux));
svmTable(2:end-2,5) = num2cell(nanstd(Yaux));
svmTable(2:end-2,6) = num2cell(round((100 - 100.*model.algosel.modelerr),1));
svmTable(2:end-2,7) = num2cell(model.algosel.svmparams(:,1));
svmTable(2:end-2,8) = num2cell(model.algosel.svmparams(:,2));              
disp(svmTable);


% ---------------------------------------------------------------------
% Storing the output data as a CSV files. This is for easier
% post-processing. All workspace data will be stored in a matlab file
% later.
writeArray2CSV = @(data,colnames,rownames,filename) writetable(array2table(data,'VariableNames',colnames,...
                                                                                'RowNames',rownames),...
                                                               filename,'WriteRowNames',true);
writeCell2CSV = @(data,colnames,rownames,filename) writetable(cell2table(data,'VariableNames',colnames,...
                                                                              'RowNames',rownames),...
                                                              filename,'WriteRowNames',true);
makeBndLabels = @(data) arrayfun(@(x) strcat('bnd_pnt_',num2str(x)),1:size(data,1),'UniformOutput',false);
colorscale  = @(data) round(255.*bsxfun(@rdivide, bsxfun(@minus, data, min(data,[],1)), range(data)));
colorscaleg = @(data) round(255.*bsxfun(@rdivide, bsxfun(@minus, data, min(data(:))), range(data(:))));

if opts.outputs.csv
    disp('-------------------------------------------------------------------------');
    disp('-> Writing the data on CSV files for post-processing.');
    % ---------------------------------------------------------------------
    writeArray2CSV(model.pbldr.Z, {'z_1','z_2'}, instlabels(subsetIndex), [rootdir 'coordinates.csv']);
    writeArray2CSV(model.sbound.Zedge, {'z_1','z_2'}, ...
                   makeBndLabels(model.sbound.Zedge), [rootdir 'bounds.csv']);
    writeArray2CSV(model.sbound.Zecorr, {'z_1','z_2'}, ...
                   makeBndLabels(model.sbound.Zecorr), [rootdir 'bounds_prunned.csv']);
    writeArray2CSV(Xraw(subsetIndex, model.featsel.idx), featlabels, instlabels(subsetIndex), [rootdir 'feature_raw.csv']);
    writeArray2CSV(X, featlabels, instlabels(subsetIndex), [rootdir 'feature_process.csv']);  
    writeArray2CSV(Yraw(subsetIndex,:), algolabels, instlabels(subsetIndex), [rootdir 'algorithm_raw.csv']);
    writeArray2CSV(Y, algolabels, instlabels(subsetIndex), [rootdir 'algorithm_process.csv']);
    writeArray2CSV(Ybin, algolabels, instlabels(subsetIndex), [rootdir 'algorithm_bin.csv']);
    writeArray2CSV(portfolio, {'Best_Algorithm'}, instlabels(subsetIndex), [rootdir 'portfolio.csv']);
    writeArray2CSV(model.algosel.Yhat, algolabels, instlabels(subsetIndex), [rootdir 'algorithm_svm.csv']);
    writeArray2CSV(model.algosel.psel, {'Best_Algorithm'}, instlabels(subsetIndex), [rootdir 'portfolio_svm.csv']);
    writeCell2CSV(model.footprint.performance(2:end,[3 5 6 8 10 11]), ...
                  model.footprint.performance(1,[3 5 6 8 10 11]),...
                  model.footprint.performance(2:end,1),...
                  [rootdir 'footprint_performance.csv']);
    writeCell2CSV(projectionMatrix(2:end,2:end), projectionMatrix(1,2:end),...
                  projectionMatrix(2:end,1), [rootdir 'projection_matrix.csv']);
    writeCell2CSV(svmTable(2:end,2:end), svmTable(1,2:end), svmTable(2:end,1), [rootdir 'svm_table.csv']);
    if opts.outputs.web
    %   writetable(array2table(parula(256), 'VariableNames', {'R','G','B'}), [rootdir 'color_table.csv']);
        writeArray2CSV(colorscale(Xraw(subsetIndex,model.featsel.idx)), featlabels, instlabels(subsetIndex), [rootdir 'feature_raw_color.csv']);
        writeArray2CSV(colorscale(Yraw(subsetIndex,:)), algolabels, instlabels(subsetIndex), [rootdir 'algorithm_raw_single_color.csv']);
        writeArray2CSV(colorscale(X), featlabels, instlabels(subsetIndex), [rootdir 'feature_process_color.csv']);
        writeArray2CSV(colorscale(Y), algolabels, instlabels(subsetIndex), [rootdir 'algorithm_process_single_color.csv']);
        writeArray2CSV(colorscaleg(Yraw(subsetIndex,:)), algolabels, instlabels(subsetIndex), [rootdir 'algorithm_raw_color.csv']);
        writeArray2CSV(colorscaleg(Y), algolabels, instlabels(subsetIndex), [rootdir 'algorithm_process_color.csv']);
    end
end
% ---------------------------------------------------------------------
% Making all the plots. First, plotting the features and performance as
% scatter plots.
if opts.outputs.png
    disp('-------------------------------------------------------------------------');
    disp('-> Producing the plots.');
    % ---------------------------------------------------------------------
    for i=1:nfeats
        clf;
        drawScatter((X(:,i)-min(X(:,i)))./range(X(:,i)), model.pbldr.Z, strrep(featlabels{i},'_',' '));
        line(model.sbound.Zedge(:,1),model.sbound.Zedge(:,2),...
                 'LineStyle', '-', ...
                 'Color', 'r');
        print(gcf,'-dpng',[rootdir 'scatter_' featlabels{i} '.png']);
    end
    % ---------------------------------------------------------------------
    Ys = log10(Yraw+1);
    Ys = (Ys-min(Ys(:)))./range(Ys(:));
    for i=1:nalgos
        clf;
        drawScatter(Ys(subsetIndex,i), model.pbldr.Z, strrep(algolabels{i},'_',' '));
        print(gcf,'-dpng',[rootdir 'scatter_' algolabels{i} '_absolute.png']);
    end
    % ---------------------------------------------------------------------
    for i=1:nalgos
        clf;
        drawScatter((Y(:,i)-min(Y(:,i)))./range(Y(:,i)), model.pbldr.Z, strrep(algolabels{i},'_',' '));
        print(gcf,'-dpng',[rootdir 'scatter_' algolabels{i} '.png']);
    end
    % ---------------------------------------------------------------------
    % Drawing the footprints for good and bad performance acording to the
    % binary measure
    for i=1:nalgos
        clf;
        drawGoodBadFootprint(model.pbldr.Z, Yfoot(:,i), model.footprint.good{i}, strrep(algolabels{i},'_',' '));
        print(gcf,'-dpng',[rootdir 'footprint_' algolabels{i} '.png']);
    end
    % ---------------------------------------------------------------------
    % Drawing the footprints as portfolio.
    clf;
    drawPortfolioFootprint(model.footprint.best, algolabels);
    print(gcf,'-dpng',[rootdir 'footprint_portfolio.png']);
    % ---------------------------------------------------------------------
    % Drawing the sources of the instances if available
    if any(issource)
        clf;
        drawSources(model.pbldr.Z, S(subsetIndex));
        print(gcf,'-dpng',[rootdir 'sources.png']);
    end
    % ---------------------------------------------------------------------
    % Drawing the SVM's predictions of good performance
    for i=1:nalgos
        clf;
        drawSVMPredictions(model.pbldr.Z, model.algosel.Yhat(:,i), strrep(algolabels{i},'_',' '));
        print(gcf,'-dpng',[rootdir 'svm_' algolabels{i} '.png']);
    end
    % ---------------------------------------------------------------------
    % Drawing the SVM's recommendations
    clf;
    drawSVMPortfolioSelections(model.pbldr.Z, model.algosel.psel, algolabels);
    print(gcf,'-dpng',[rootdir 'svm_portfolio.png']);
end
% -------------------------------------------------------------------------
model.opts = opts;
disp('-------------------------------------------------------------------------');
disp('-> Storing the raw MATLAB results for post-processing and/or debugging.');
save([rootdir 'model.mat'],'-struct','model'); % Save the main results
save([rootdir 'workspace.mat']); % Save the full workspace for debugging
disp(['-> Completed! Elapsed time: ' num2str(toc(startProcess)) 's']);
disp('EOF:SUCCESS');
end
