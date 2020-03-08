% -------------------------------------------------------------------------
% csvscript.m
% -------------------------------------------------------------------------
%
% By: Mario Andres Munoz Acosta
%     School of Mathematics and Statistics
%     The University of Melbourne
%     Australia
%     2020
%
% -------------------------------------------------------------------------

disp('=========================================================================');
disp('-> Writing the data on CSV files for posterior analysis.');
% -------------------------------------------------------------------------
for i=1:nalgos
    if isfield(model.trace.best{i},'polygon') && ~isempty(model.trace.best{i}.polygon)
        writeArray2CSV(model.trace.best{i}.polygon.Vertices, {'z_1','z_2'},...
                       makeBndLabels(model.trace.best{i}.polygon.Vertices),...
                       [rootdir 'footprint_' algolabels{i} '_best.csv']);
    end
    if isfield(model.trace.good{i},'polygon') && ~isempty(model.trace.good{i}.polygon)
        writeArray2CSV(model.trace.good{i}.polygon.Vertices, {'z_1','z_2'},...
                       makeBndLabels(model.trace.good{i}.polygon.Vertices),...
                       [rootdir 'footprint_' algolabels{i} '_good.csv']);
    end
    if isfield(model.trace.bad{i},'polygon') && ~isempty(model.trace.bad{i}.polygon)
        writeArray2CSV(model.trace.bad{i}.polygon.Vertices, {'z_1','z_2'},...
                       makeBndLabels(model.trace.bad{i}.polygon.Vertices),...
                       [rootdir 'footprint_' algolabels{i} '_bad.csv']);
    end
end

writeArray2CSV(model.pilot.Z, {'z_1','z_2'}, instlabels(subsetIndex), [rootdir 'coordinates.csv']);
writeArray2CSV(model.cloist.Zedge, {'z_1','z_2'}, makeBndLabels(model.cloist.Zedge), [rootdir 'bounds.csv']);
writeArray2CSV(model.cloist.Zecorr, {'z_1','z_2'}, makeBndLabels(model.cloist.Zecorr), [rootdir 'bounds_prunned.csv']);
writeArray2CSV(Xraw(subsetIndex, model.featsel.idx), featlabels, instlabels(subsetIndex), [rootdir 'feature_raw.csv']);
writeArray2CSV(X, featlabels, instlabels(subsetIndex), [rootdir 'feature_process.csv']);
writeArray2CSV(Yraw(subsetIndex,:), algolabels, instlabels(subsetIndex), [rootdir 'algorithm_raw.csv']);
writeArray2CSV(Y, algolabels, instlabels(subsetIndex), [rootdir 'algorithm_process.csv']);
writeArray2CSV(Ybin, algolabels, instlabels(subsetIndex), [rootdir 'algorithm_bin.csv']);
writeArray2CSV(numGoodAlgos(subsetIndex), {'NumGoodAlgos'}, instlabels(subsetIndex), [rootdir 'good_algos.csv']);
writeArray2CSV(beta, {'IsBetaEasy'}, instlabels(subsetIndex), [rootdir 'beta_easy.csv']);
writeArray2CSV(P, {'Best_Algorithm'}, instlabels(subsetIndex), [rootdir 'portfolio.csv']);
writeArray2CSV(model.pythia.Yhat, algolabels, instlabels(subsetIndex), [rootdir 'algorithm_svm.csv']);
writeArray2CSV(model.pythia.selection0, {'Best_Algorithm'}, instlabels(subsetIndex), [rootdir 'portfolio_svm.csv']);
writeCell2CSV(model.trace.summary(2:end,[3 5 6 8 10 11]), ...
              model.trace.summary(1,[3 5 6 8 10 11]),...
              model.trace.summary(2:end,1),...
              [rootdir 'footprint_performance.csv']);
writeCell2CSV(model.pilot.summary(2:end,2:end), model.pilot.summary(1,2:end),...
              model.pilot.summary(2:end,1), [rootdir 'projection_matrix.csv']);
writeCell2CSV(model.pythia.summary(2:end,2:end), ...
              model.pythia.summary(1,2:end), ...
              model.pythia.summary(2:end,1), [rootdir 'svm_table.csv']);
