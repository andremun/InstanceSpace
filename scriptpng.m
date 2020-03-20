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
nfeats = size(model.data.X,2);
Xaux = (model.data.X-min(model.data.X,[],1))./range(model.data.X,1);
Yind = (model.data.Y-min(model.data.Y,[],1))./range(model.data.Y,1);
Yglb = log10(model.data.Yraw+1);
Yglb = (Yglb-min(Yglb(:)))./range(Yglb(:));
if opts.trace.usesim
    Yfoot = model.pythia.Yhat;
    Pfoot = model.pythia.selection0;
else
    Yfoot = model.data.Ybin;
    Pfoot = model.data.P;
end
% -------------------------------------------------------------------------
disp('=========================================================================');
disp('-> Producing the plots.');
% -------------------------------------------------------------------------
% Drawing feature plots
for i=1:nfeats
    clf;
    drawScatter(model.pilot.Z, Xaux(:,i),...
                strrep(model.data.featlabels{i},'_',' '));
    % line(model.cloist.Zedge(:,1), model.cloist.Zedge(:,2), 'LineStyle', '-', 'Color', 'r');
    print(gcf,'-dpng',[rootdir 'distribution_feature_' model.data.featlabels{i} '.png']);
end
% -------------------------------------------------------------------------
% Drawing algorithm performance/footprint plots
for i=1:nalgos
    % Actual performance, normalized globaly
    clf;
    drawScatter(model.pilot.Z, Yglb(:,i), ...
                strrep(model.data.algolabels{i},'_',' '));
    print(gcf,'-dpng',[rootdir 'distribution_performance_global_normalized_' model.data.algolabels{i} '.png']);
    % Actual performance, normalized individualy
    clf;
    drawScatter(model.pilot.Z, Yind(:,i), ...
                strrep(model.data.algolabels{i},'_',' '));
    print(gcf,'-dpng',[rootdir 'distribution_performance_individual_normalized_' model.data.algolabels{i} '.png']);
    % Actual binary performance
    clf;
    drawBinaryPerformance(model.pilot.Z, model.data.Ybin(:,i), ...
                          strrep(model.data.algolabels{i},'_',' '));
    print(gcf,'-dpng',[rootdir 'binary_performance_' model.data.algolabels{i} '.png']);
    % Drawing the SVM's predictions of good performance
    clf;
    drawBinaryPerformance(model.pilot.Z, model.pythia.Yhat(:,i), ...
                          strrep(model.data.algolabels{i},'_',' '));
    print(gcf,'-dpng',[rootdir 'binary_svm_' model.data.algolabels{i} '.png']);
    % Drawing the footprints for good and bad performance acording to the
    % binary measure 
    clf;
    drawGoodBadFootprint(model.pilot.Z, ...
                         model.trace.good{i}, ...
                         model.trace.bad{i}, ...
                         Yfoot(:,i), ...
                         strrep(model.data.algolabels{i},'_',' '));
    print(gcf,'-dpng',[rootdir 'footprint_' model.data.algolabels{i} '.png']);
end
% ---------------------------------------------------------------------
% Plotting the number of good algos
clf;
drawScatter(model.pilot.Z, model.data.numGoodAlgos./nalgos, 'Percentage of good algorithms');
print(gcf,'-dpng',[rootdir 'distribution_number_good_algos.png']);
% ---------------------------------------------------------------------
% Drawing the algorithm performance
clf;
drawPortfolioSelections(model.pilot.Z, model.data.P, model.data.algolabels, 'Best algorithm');
print(gcf,'-dpng',[rootdir 'distribution_portfolio.png']);
% ---------------------------------------------------------------------
% Drawing the SVM's recommendations
clf;
drawPortfolioSelections(model.pilot.Z, model.pythia.selection0, model.data.algolabels, 'Predicted best algorithm');
print(gcf,'-dpng',[rootdir 'distribution_svm_portfolio.png']);
% ---------------------------------------------------------------------
% Drawing the footprints as portfolio.
clf;
drawPortfolioFootprint(model.pilot.Z, model.trace.best, Pfoot, model.data.algolabels);
print(gcf,'-dpng',[rootdir 'footprint_portfolio.png']);
% ---------------------------------------------------------------------
% Plotting the model.data.beta score
clf;
drawBinaryPerformance(model.pilot.Z, model.data.beta, '\beta score');
print(gcf,'-dpng',[rootdir 'distribution_beta_score.png']);
% ---------------------------------------------------------------------
% Drawing the sources of the instances if available
if any(issource)
    clf;
    drawSources(model.pilot.Z, model.data.S);
    print(gcf,'-dpng',[rootdir 'distribution_sources.png']);
end