function scriptpng(container,rootdir)
% -------------------------------------------------------------------------
% pgnscript.m
% -------------------------------------------------------------------------
%
% By: Mario Andres Munoz Acosta
%     School of Mathematics and Statistics
%     The University of Melbourne
%     Australia
%     2020
%
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Preliminaries
scriptfcn;
colormap('parula');
nfeats = size(container.data.X,2);
nalgos = size(container.data.Y,2);
Xaux = (container.data.X-min(container.data.X,[],1))./range(container.data.X,1);
Yind = (container.data.Yraw-min(container.data.Yraw,[],1))./range(container.data.Yraw,1);
Yglb = log10(container.data.Yraw+1);
Yglb = (Yglb-min(Yglb(:)))./range(Yglb(:));
if container.opts.trace.usesim
    Yfoot = container.pythia.Yhat;
    Pfoot = container.pythia.selection0;
else
    Yfoot = container.data.Ybin;
    Pfoot = container.data.P;
end
% -------------------------------------------------------------------------
disp('=========================================================================');
disp('-> Producing the plots.');
% -------------------------------------------------------------------------
% Drawing feature plots
for i=1:nfeats
    clf;
    drawScatter(container.pilot.Z, Xaux(:,i),...
                strrep(container.data.featlabels{i},'_',' '));
    % line(model.cloist.Zedge(:,1), model.cloist.Zedge(:,2), 'LineStyle', '-', 'Color', 'r');
    print(gcf,'-dpng',[rootdir 'distribution_feature_' container.data.featlabels{i} '.png']);
end
% -------------------------------------------------------------------------
% Drawing algorithm performance/footprint plots
for i=1:nalgos
    % Actual performance, normalized globaly
    clf;
    drawScatter(container.pilot.Z, Yglb(:,i), ...
                strrep(container.data.algolabels{i},'_',' '));
    print(gcf,'-dpng',[rootdir 'distribution_performance_global_normalized_' container.data.algolabels{i} '.png']);
    % Actual performance, normalized individualy
    clf;
    drawScatter(container.pilot.Z, Yind(:,i), ...
                strrep(container.data.algolabels{i},'_',' '));
    print(gcf,'-dpng',[rootdir 'distribution_performance_individual_normalized_' container.data.algolabels{i} '.png']);
    % Actual binary performance
    clf;
    drawBinaryPerformance(container.pilot.Z, container.data.Ybin(:,i), ...
                          strrep(container.data.algolabels{i},'_',' '));
    print(gcf,'-dpng',[rootdir 'binary_performance_' container.data.algolabels{i} '.png']);
    % Drawing the SVM's predictions of good performance
    try
        clf;
        drawBinaryPerformance(container.pilot.Z, container.pythia.Yhat(:,i), ...
                              strrep(container.data.algolabels{i},'_',' '));
        print(gcf,'-dpng',[rootdir 'binary_svm_' container.data.algolabels{i} '.png']);
    catch
        disp('No SVM model has been trained');
    end
    % Drawing the footprints for good and bad performance acording to the
    % binary measure 
    try 
        clf;
        drawGoodBadFootprint(container.pilot.Z, ...
                             container.trace.good{i}, ...
                             Yfoot(:,i), ...
                             strrep(container.data.algolabels{i},'_',' '));
        print(gcf,'-dpng',[rootdir 'footprint_' container.data.algolabels{i} '.png']);
    catch
        disp('No Footprint has been calculated');
    end
end
% ---------------------------------------------------------------------
% Plotting the number of good algos
clf;
drawScatter(container.pilot.Z, container.data.numGoodAlgos./nalgos, 'Percentage of good algorithms');
print(gcf,'-dpng',[rootdir 'distribution_number_good_algos.png']);
% ---------------------------------------------------------------------
% Drawing the algorithm performance
clf;
drawPortfolioSelections(container.pilot.Z, container.data.P, container.data.algolabels, 'Best algorithm');
print(gcf,'-dpng',[rootdir 'distribution_portfolio.png']);
% ---------------------------------------------------------------------
% Drawing the SVM's recommendations
clf;
drawPortfolioSelections(container.pilot.Z, container.pythia.selection0, container.data.algolabels, 'Predicted best algorithm');
print(gcf,'-dpng',[rootdir 'distribution_svm_portfolio.png']);
% ---------------------------------------------------------------------
% Drawing the footprints as portfolio.
clf;
drawPortfolioFootprint(container.pilot.Z, container.trace.best, Pfoot, container.data.algolabels);
print(gcf,'-dpng',[rootdir 'footprint_portfolio.png']);
% ---------------------------------------------------------------------
% Plotting the model.data.beta score
clf;
drawBinaryPerformance(container.pilot.Z, container.data.beta, '\beta score');
print(gcf,'-dpng',[rootdir 'distribution_beta_score.png']);
% ---------------------------------------------------------------------
% Drawing the sources of the instances if available
if isfield(container.data,'S')
    clf;
    drawSources(container.pilot.Z, container.data.S);
    print(gcf,'-dpng',[rootdir 'distribution_sources.png']);
end