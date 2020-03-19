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

nfeats = size(X,2);
disp('=========================================================================');
disp('-> Producing the plots.');
% -------------------------------------------------------------------------
% Drawing feature plots
for i=1:nfeats
    clf;
    drawScatter(model.pilot.Z, (X(:,i)-min(X(:,i)))./range(X(:,i)), strrep(featlabels{i},'_',' '));
    % line(model.cloist.Zedge(:,1), model.cloist.Zedge(:,2), 'LineStyle', '-', 'Color', 'r');
    print(gcf,'-dpng',[rootdir 'distribution_feature_' featlabels{i} '.png']);
end
% -------------------------------------------------------------------------
% Drawing algorithm performance/footprint plots
Ys = log10(Yraw+1);
Ys = (Ys-min(Ys(:)))./range(Ys(:));
if opts.trace.usesim
    Yfoot = model.pythia.Yhat;
    Pfoot = model.pythia.selection0;
else
    Yfoot = Ybin;
    Pfoot = P;
end
for i=1:nalgos
    % Actual performance, normalized globaly
    clf;
    drawScatter(model.pilot.Z, Ys(:,i), strrep(algolabels{i},'_',' '));
    print(gcf,'-dpng',[rootdir 'distribution_performance_global_normalized_' algolabels{i} '.png']);
    % Actual performance, normalized individualy
    clf;
    drawScatter(model.pilot.Z, (Y(:,i)-min(Y(:,i)))./range(Y(:,i)), strrep(algolabels{i},'_',' '));
    print(gcf,'-dpng',[rootdir 'distribution_performance_individual_normalized_' algolabels{i} '.png']);
    % Actual binary performance
    clf;
    drawBinaryPerformance(model.pilot.Z, Ybin(:,i), strrep(algolabels{i},'_',' '));
    print(gcf,'-dpng',[rootdir 'binary_performance_' algolabels{i} '.png']);
    % Drawing the SVM's predictions of good performance
    clf;
    drawBinaryPerformance(model.pilot.Z, model.pythia.Yhat(:,i), strrep(algolabels{i},'_',' '));
    print(gcf,'-dpng',[rootdir 'binary_svm_' algolabels{i} '.png']);
    % Drawing the footprints for good and bad performance acording to the
    % binary measure 
    clf;
    drawGoodBadFootprint(model.pilot.Z, model.trace.good{i}, model.trace.bad{i}, Yfoot(:,i), strrep(algolabels{i},'_',' '));
    print(gcf,'-dpng',[rootdir 'footprint_' algolabels{i} '.png']);
end
% ---------------------------------------------------------------------
% Plotting the number of good algos
clf;
drawScatter(model.pilot.Z, numGoodAlgos./nalgos, 'Percentage of good algorithms');
print(gcf,'-dpng',[rootdir 'distribution_number_good_algos.png']);
% ---------------------------------------------------------------------
% Drawing the algorithm performance
clf;
drawPortfolioSelections(model.pilot.Z, P, algolabels, 'Best algorithm');
print(gcf,'-dpng',[rootdir 'distribution_portfolio.png']);
% ---------------------------------------------------------------------
% Drawing the SVM's recommendations
clf;
drawPortfolioSelections(model.pilot.Z, model.pythia.selection0, algolabels, 'Predicted best algorithm');
print(gcf,'-dpng',[rootdir 'distribution_svm_portfolio.png']);
% ---------------------------------------------------------------------
% Drawing the footprints as portfolio.
clf;
drawPortfolioFootprint(model.pilot.Z, model.trace.best, Pfoot, algolabels);
print(gcf,'-dpng',[rootdir 'footprint_portfolio.png']);
% ---------------------------------------------------------------------
% Plotting the beta score
clf;
drawBinaryPerformance(model.pilot.Z, beta, '\beta score');
print(gcf,'-dpng',[rootdir 'distribution_beta_score.png']);
% ---------------------------------------------------------------------
% Drawing the sources of the instances if available
if any(issource)
    clf;
    drawSources(model.pilot.Z, S);
    print(gcf,'-dpng',[rootdir 'distribution_sources.png']);
end