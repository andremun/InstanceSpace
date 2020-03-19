% -------------------------------------------------------------------------
% webscript.m
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
disp('-> Writing the data for the web interfase.');
% -------------------------------------------------------------------------
writetable(array2table(parula(256), 'VariableNames', {'R','G','B'}), [rootdir 'color_table.csv']);
writeArray2CSV(colorscale(Xraw(:,model.featsel.idx)), featlabels, instlabels, [rootdir 'feature_raw_color.csv']);
writeArray2CSV(colorscale(Yraw), algolabels, instlabels, [rootdir 'algorithm_raw_single_color.csv']);
writeArray2CSV(colorscale(X), featlabels, instlabels, [rootdir 'feature_process_color.csv']);
writeArray2CSV(colorscale(Y), algolabels, instlabels, [rootdir 'algorithm_process_single_color.csv']);
writeArray2CSV(colorscaleg(Yraw), algolabels, instlabels, [rootdir 'algorithm_raw_color.csv']);
writeArray2CSV(colorscaleg(Y), algolabels, instlabels, [rootdir 'algorithm_process_color.csv']);
writeArray2CSV(colorscaleg(numGoodAlgos),  {'NumGoodAlgos'}, instlabels, [rootdir 'good_algos_color.csv']);