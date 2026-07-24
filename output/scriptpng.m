% -------------------------------------------------------------------------
% Instance Space Analysis (ISA) Toolkit
% Copyright (c) 2026 Mario Andres Munoz Acosta and contributors
% School of Computing and Information Systems
% The University of Melbourne, Australia
%
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
% License: https://polyformproject.org/licenses/noncommercial/1.0.0/
%
% You may use, modify, and redistribute this software for non-commercial
% research and educational purposes only. Commercial use requires prior
% written permission. See the LICENSE file for full terms.
%
% Reference:
%   Simpson, C., Munoz, M.A., Kandanaarachchi, S. & Campello, R.J.G.B.
%   (2025). ISA3: A 3-dimensional expansion of Instance Space Analysis.
%   Machine Learning, 114, 240. https://doi.org/10.1007/s10994-025-06871-5
%
%   Smith-Miles, K. & Munoz, M.A. (2023). Instance Space Analysis for
%   Algorithm Testing. ACM Computing Surveys, 55(12), Article 255.
%   https://doi.org/10.1145/3572895
% -------------------------------------------------------------------------
function scriptpng(container,rootdir)
% scriptpng  Write a model or explore() result to PNG figures in rootdir.
%
%   scriptpng(container,rootdir)
%
%   container - model struct from buildIS/InstanceSpace.build(), or a
%               testResults entry from exploreIS/InstanceSpace.explore()
%   rootdir   - destination directory (trailing slash required)
%
%   Produces per-feature and per-algorithm distribution plots, portfolio
%   selection and footprint plots, using scriptfcn.m's drawing helpers.
%   Renders in 3D and applies the optimised camera viewpoint
%   (container.pilot.viewpoint, see PILOTviewpoint.m) when the projection
%   is 3D and one was computed.

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
Yfoot = container.data.Ybin;
Pfoot = container.data.P;
% -------------------------------------------------------------------------
% Optimised 3D camera viewpoint(s) (spec §5.2), if computed. Absent for 2D
% projections; resolveViewAngle then returns [] and applyView is a no-op
% in 2D (it only calls the default view(3) as a fallback for a 3D plot
% that has no computed viewpoint to use).
if isfield(container.pilot,'viewpoint')
    viewpoint = container.pilot.viewpoint;
else
    viewpoint = [];
end
globalView = resolveViewAngle(viewpoint, []); % feature/portfolio-level plots
% -------------------------------------------------------------------------
fprintf('[OUTPUT] Producing the plots.\n');
% -------------------------------------------------------------------------
% Drawing feature plots
for i=1:nfeats
    clf;
    drawScatter(container.pilot.Z, Xaux(:,i),...
                strrep(container.data.featlabels{i},'_',' '), globalView);
    % line(model.cloist.Zedge(:,1), model.cloist.Zedge(:,2), 'LineStyle', '-', 'Color', 'r');
    print(gcf,'-dpng',[rootdir 'distribution_feature_' container.data.featlabels{i} '.png']);
end
% -------------------------------------------------------------------------
% Drawing algorithm performance/footprint plots
for i=1:nalgos
    algoView = resolveViewAngle(viewpoint, i);
    % Actual performance, normalized globaly
    clf;
    drawScatter(container.pilot.Z, Yglb(:,i), ...
                strrep(container.data.algolabels{i},'_',' '), algoView);
    print(gcf,'-dpng',[rootdir 'distribution_performance_global_normalized_' container.data.algolabels{i} '.png']);
    % Actual performance, normalized individualy
    clf;
    drawScatter(container.pilot.Z, Yind(:,i), ...
                strrep(container.data.algolabels{i},'_',' '), algoView);
    print(gcf,'-dpng',[rootdir 'distribution_performance_individual_normalized_' container.data.algolabels{i} '.png']);
    % Actual binary performance
    try
        clf;
        drawBinaryPerformance(container.pilot.Z, container.data.Ybin(:,i), ...
                              strrep(container.data.algolabels{i},'_',' '), algoView);
        print(gcf,'-dpng',[rootdir 'binary_performance_' container.data.algolabels{i} '.png']);
    catch
        fprintf('[OUTPUT] No binary performance has been calculated.\n');
    end
    % Drawing the classifier's predictions of good performance
    try
        clf;
        drawBinaryPerformance(container.pilot.Z, container.pythia.Yhat(:,i), ...
                              strrep(container.data.algolabels{i},'_',' '), algoView);
        print(gcf,'-dpng',[rootdir 'binary_classifier_' container.data.algolabels{i} '.png']);
    catch
        fprintf('[OUTPUT] No classifier predictions are available.\n');
    end
    % Drawing the footprints for good and bad performance according to the
    % binary measure
    try
        clf;
        drawGoodBadFootprint(container.pilot.Z, ...
                             container.trace.good{i}, ...
                             Yfoot(:,i), ...
                             strrep(container.data.algolabels{i},'_',' '), algoView);
        print(gcf,'-dpng',[rootdir 'footprint_' container.data.algolabels{i} '.png']);
    catch
        fprintf('[OUTPUT] No Footprint has been calculated.\n');
    end
end
% ---------------------------------------------------------------------
% Plotting the number of good algos
clf;
drawScatter(container.pilot.Z, container.data.numGoodAlgos./nalgos, 'Percentage of good algorithms', globalView);
print(gcf,'-dpng',[rootdir 'distribution_number_good_algos.png']);
% ---------------------------------------------------------------------
% Drawing the algorithm performance
clf;
drawPortfolioSelections(container.pilot.Z, container.data.P, container.data.algolabels, 'Best algorithm', globalView);
print(gcf,'-dpng',[rootdir 'distribution_portfolio.png']);
% ---------------------------------------------------------------------
% Drawing the SVM's recommendations
clf;
drawPortfolioSelections(container.pilot.Z, container.pythia.selection0, container.data.algolabels, 'Predicted best algorithm', globalView);
print(gcf,'-dpng',[rootdir 'distribution_svm_portfolio.png']);
% ---------------------------------------------------------------------
% Drawing the footprints as portfolio.
clf;
drawPortfolioFootprint(container.pilot.Z, container.trace.best, Pfoot, container.data.algolabels, globalView);
print(gcf,'-dpng',[rootdir 'footprint_portfolio.png']);
% ---------------------------------------------------------------------
% Plotting the model.data.beta score
clf;
drawBinaryPerformance(container.pilot.Z, container.data.beta, '\beta score', globalView);
print(gcf,'-dpng',[rootdir 'distribution_beta_score.png']);
% ---------------------------------------------------------------------
% Drawing the sources of the instances if available
if isfield(container.data,'S')
    clf;
    drawSources(container.pilot.Z, container.data.S, globalView);
    print(gcf,'-dpng',[rootdir 'distribution_sources.png']);
end