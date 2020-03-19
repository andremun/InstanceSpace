function out = testIS(rootdir)
% -------------------------------------------------------------------------
% testIS.m
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
scriptdisc('testIS.m');
% -------------------------------------------------------------------------
% Collect all the data from the files
disp(['Root Directory: ' rootdir]);
modelfile = [rootdir 'model.mat'];
datafile = [rootdir 'metadata.csv'];
optsfile = [rootdir 'options.json'];
if ~isfile(modelfile) || ~isfile(datafile) || ~isfile(optsfile)
    error(['Please place the datafiles in the directory ''' rootdir '''']);
end
opts = jsondecode(fileread(optsfile));
disp('-------------------------------------------------------------------------');
disp('Listing options used:');
optfields = fieldnames(opts);
for i = 1:length(optfields)
    disp(optfields{i});
    disp(opts.(optfields{i}));
end
disp('-------------------------------------------------------------------------');
disp('-> Loading the data');
model = load(modelfile);
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
nalgos = size(Y,2);
ninst = size(X,1);
% -------------------------------------------------------------------------
% Storing the raw data for further processing, e.g., graphs
Xraw = X;
Yraw = Y;
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
MaxPerf = false;
if isfield(model.opts.perf, 'MaxPerf')
    MaxPerf = model.opts.perf.MaxPerf;
elseif  isfield(model.opts.perf, 'MaxMin')
    MaxPerf = model.opts.perf.MaxMin;
else
    warning('Can not find parameter "MaxPerf" in the trained model. We are assuming that performance metric is needed to be minimized.');
end
if MaxPerf
    Yaux = Y;
    Yaux(isnan(Yaux)) = -Inf;
    [rankPerf,rankAlgo] = sort(Yaux,2,'descend');
    bestPerformace = rankPerf(:,1);
    P = rankAlgo(:,1);
    if model.opts.perf.AbsPerf
        Ybin = Y>=model.opts.perf.epsilon;
        msg = [msg 'higher than ' num2str(model.opts.perf.epsilon)];
    else
        Ybin = bsxfun(@ge,Y,(1-model.opts.perf.epsilon).*bestPerformace); % One is good, zero is bad
        msg = [msg 'within ' num2str(round(100.*model.opts.perf.epsilon)) '% of the best.'];
    end
else
    Yaux = Y;
    Yaux(isnan(Yaux)) = Inf;
    [rankPerf,rankAlgo] = sort(Yaux,2,'ascend');
    bestPerformace = rankPerf(:,1);
    P = rankAlgo(:,1);
    if model.opts.perf.AbsPerf
        Ybin = Y<=model.opts.perf.epsilon;
        msg = [msg 'less than ' num2str(model.opts.perf.epsilon)];
    else
        Ybin = bsxfun(@le,Y,(1+model.opts.perf.epsilon).*bestPerformace);
        msg = [msg 'within ' num2str(round(100.*model.opts.perf.epsilon)) '% of the best.'];
    end
end
disp(msg);
numGoodAlgos = sum(Ybin,2);
beta = numGoodAlgos>model.opts.general.betaThreshold*nalgos;

% ---------------------------------------------------------------------
% Automated pre-processing
if model.opts.auto.preproc && model.opts.bound.flag
    disp('-------------------------------------------------------------------------');
    disp('-> Auto-pre-processing. Bounding outliers, scaling and normalizing the data.');
    % Eliminate extreme outliers, i.e., any point that exceedes 5 times the
    % inter quantile range, by bounding them to that value.
    disp('-> Removing extreme outliers from the feature values.');
    himask = bsxfun(@gt,X,model.bound.hibound);
    lomask = bsxfun(@lt,X,model.bound.lobound);
    X = X.*~(himask | lomask) + bsxfun(@times,himask,model.bound.hibound) + ...
                                bsxfun(@times,lomask,model.bound.lobound);
end

if model.opts.auto.preproc && model.opts.norm.flag
    % Normalize the data using Box-Cox and out.pilot.Z-transformations
    disp('-> Auto-normalizing the data.');
    X = bsxfun(@minus,X,model.norm.minX)+1;
    X = bsxfun(@rdivide,bsxfun(@power,X,model.norm.lambdaX)-1,model.norm.lambdaX);
    X = bsxfun(@rdivide,bsxfun(@minus,X,model.norm.muX),model.norm.sigmaX);
    
    Y(Y==0) = eps; % Assumes that Y is always positive and higher than 1e-16
    Y = bsxfun(@rdivide,bsxfun(@power,Y,model.norm.lambdaY)-1,model.norm.lambdaY);
    Y = bsxfun(@rdivide,bsxfun(@minus,Y,model.norm.muY),model.norm.sigmaY);
end
% ---------------------------------------------------------------------
% This is the final subset of features.
X = X(:,model.featsel.idx);
nfeats = size(X,2);
featlabels = strrep(varlabels(isfeat),'feature_','');
featlabels = featlabels(model.featsel.idx);
algolabels = strrep(varlabels(isalgo),'algo_','');
% ---------------------------------------------------------------------
%  Calculate the two dimensional projection using the PBLDR algorithm
%  (Munoz et al. Mach Learn 2018)
out.pilot.Z = X*model.pilot.A';
% -------------------------------------------------------------------------
% Algorithm selection. Fit a model that would separate the space into
% classes of good and bad performance. 
out.pythia.Yhat = false(size(Ybin));
out.pythia.Pr0hat = 0.*Ybin;
cvcmat = zeros(nalgos,4);
for i=1:nalgos
    [out.pythia.Yhat(:,i),aux] = model.pythia.svm{i}.predict(out.pilot.Z);
    out.pythia.Pr0hat(:,i) = aux(:,1);
    aux = confusionmat(Ybin(:,i),out.pythia.Yhat(:,i));
    cvcmat(i,:) = aux(:);
end
tn = cvcmat(:,1);
fp = cvcmat(:,3);
fn = cvcmat(:,2);
tp = cvcmat(:,4);
precision = tp./(tp+fp);
recall = tp./(tp+fn);
accuracy = (tp+tn)./sum(cvcmat(1,:));

[mostprecise,psel] = max(bsxfun(@times,out.pythia.Yhat,model.pythia.precision'),[],2);
pselfull = psel;
psel(mostprecise<=0) = 0;
[~,betterdefault] = max(mean(Ybin));
pselfull(mostprecise<=0) = betterdefault;

svmselections = bsxfun(@eq,psel,1:nalgos);
selselections = bsxfun(@eq,pselfull,1:nalgos);
avgperf = mean(Yraw);
stdperf = std(Yraw);
Yselector = Yraw;
Yfull = Yraw;
Ysvms = Yraw;
Yselector(~svmselections) = NaN;
Yfull(~selselections) = NaN;
Ysvms(~out.pythia.Yhat) = NaN;

pgood = mean(any( Ybin & selselections,2));
fb = sum(any( Ybin & ~svmselections,2));
fg = sum(any(~Ybin &  svmselections,2));
tg = sum(any( Ybin &  svmselections,2));
precisionsel = tg./(tg+fg);
recallsel = tg./(tg+fb);

out.pythia.summary = cell(nalgos+3, 9);
out.pythia.summary{1,1} = 'Algorithms ';
out.pythia.summary(2:end-2, 1) = algolabels;
out.pythia.summary(end-1:end, 1) = {'Oracle','Selector'};
out.pythia.summary(1, 2:9) = {'Avg_Perf_all_instances';
                              'Std_Perf_all_instances';
                              'Probability_of_good';
                              'Avg_Perf_selected_instances';
                              'Std_Perf_selected_instances';
                              'model_accuracy';
                              'model_precision';
                              'model_recall'};
out.pythia.summary(2:end, 2) = num2cell(round([avgperf nanmean(bestPerformace) nanmean(Yfull(:))],3));
out.pythia.summary(2:end, 3) = num2cell(round([stdperf nanstd(bestPerformace) nanstd(Yfull(:))],3));
out.pythia.summary(2:end, 4) = num2cell(round([mean(Ybin) 1 pgood],3));
out.pythia.summary(2:end, 5) = num2cell(round([nanmean(Ysvms) NaN nanmean(Yselector(:))],3));
out.pythia.summary(2:end, 6) = num2cell(round([nanstd(Ysvms) NaN nanstd(Yselector(:))],3));
out.pythia.summary(2:end, 7) = num2cell(round(100.*[accuracy' NaN NaN],1));
out.pythia.summary(2:end, 8) = num2cell(round(100.*[precision' NaN precisionsel],1));
out.pythia.summary(2:end, 9) = num2cell(round(100.*[recall' NaN recallsel],1));
out.pythia.summary(cellfun(@(x) all(isnan(x)),out.pythia.summary)) = {[]}; % Clean up. Not really needed
disp('-> Completed! Performance of the models:');
disp(' ');
disp(out.pythia.summary);

% -------------------------------------------------------------------------
% Writing the results
scriptfcn;

if opts.outputs.csv
    disp('-------------------------------------------------------------------------');
    disp('-> Writing the data on CSV files for post-processing.');
    % ---------------------------------------------------------------------
    writeArray2CSV(out.pilot.Z, {'z_1','z_2'}, instlabels, [rootdir 'coordinates_test.csv']);
    writeArray2CSV(Xraw(:, model.featsel.idx), featlabels, instlabels, [rootdir 'feature_test_raw.csv']);
    writeArray2CSV(X, featlabels, instlabels, [rootdir 'feature_test_process.csv']);  
    writeArray2CSV(Yraw, algolabels, instlabels, [rootdir 'algorithm_test_raw.csv']);
    writeArray2CSV(Y, algolabels, instlabels, [rootdir 'algorithm_test_process.csv']);
    writeArray2CSV(Ybin, algolabels, instlabels, [rootdir 'algorithm_test_bin.csv']);
    writeArray2CSV(numGoodAlgos, {'NumGoodAlgos'}, instlabels, [rootdir 'good_algos_test.csv']);
    writeArray2CSV(beta, {'IsBetaEasy'}, instlabels, [rootdir 'beta_easy_test.csv']);
    writeArray2CSV(P, {'Best_Algorithm'}, instlabels, [rootdir 'portfolio_test.csv']);
    writeArray2CSV(out.pythia.Yhat, algolabels, instlabels, [rootdir 'algorithm_test_svm.csv']);
    writeArray2CSV(psel, {'Best_Algorithm'}, instlabels, [rootdir 'portfolio_test_svm.csv']);
    writeCell2CSV(out.pythia.summary(2:end,2:end), out.pythia.summary(1,2:end), out.pythia.summary(2:end,1), [rootdir 'svm_test_table.csv']);
    if opts.outputs.web
    %   writetable(array2table(parula(256), 'VariableNames', {'R','G','B'}), [rootdir 'color_table.csv']);
        writeArray2CSV(colorscale(Xraw(:,model.featsel.idx)), featlabels, instlabels, [rootdir 'feature_test_raw_color.csv']);
        writeArray2CSV(colorscale(Yraw), algolabels, instlabels, [rootdir 'algorithm_test_raw_single_color.csv']);
        writeArray2CSV(colorscale(X), featlabels, instlabels, [rootdir 'feature_test_process_color.csv']);
        writeArray2CSV(colorscale(Y), algolabels, instlabels, [rootdir 'algorithm_test_process_single_color.csv']);
        writeArray2CSV(colorscaleg(Yraw), algolabels, instlabels, [rootdir 'algorithm_test_raw_color.csv']);
        writeArray2CSV(colorscaleg(Y), algolabels, instlabels, [rootdir 'algorithm_test_process_color.csv']);
        writeArray2CSV(colorscael(numGoodAlgos), {'NumGoodAlgos'}, instlabels, [rootdir 'good_algos_test_color.csv']);
    end
end

if opts.outputs.png
    disp('-------------------------------------------------------------------------');
    disp('-> Producing the plots.');
    % ---------------------------------------------------------------------
    for i=1:nfeats
        clf;
        drawScatter(out.pilot.Z, (X(:,i)-min(X(:,i)))./range(X(:,i)), strrep(featlabels{i},'_',' '));
        % line(model.sbound.out.pilot.Zedge(:,1),model.sbound.out.pilot.Zedge(:,2),'LineStyle', '-', 'Color', 'r');
        print(gcf,'-dpng',[rootdir 'distribution_test_feature_' featlabels{i} '.png']);
    end
    % ---------------------------------------------------------------------
    Ys = log10(Yraw+1);
    Ys = (Ys-min(Ys(:)))./range(Ys(:));
    for i=1:nalgos
        clf;
        drawScatter(out.pilot.Z, Ys(:,i), strrep(algolabels{i},'_',' '));
        print(gcf,'-dpng',[rootdir 'distribution_test_performance_global_normalized_' algolabels{i} '.png']);
    end
    % ---------------------------------------------------------------------
    for i=1:nalgos
        clf;
        drawScatter(out.pilot.Z, (Y(:,i)-min(Y(:,i)))./range(Y(:,i)), strrep(algolabels{i},'_',' '));
        print(gcf,'-dpng',[rootdir 'distribution_test_performance_individual_normalized_' algolabels{i} '.png']);
    end
    % ---------------------------------------------------------------------
    % Drawing the sources of the instances if available
    if any(issource)
        clf;
        drawSources(out.pilot.Z, S);
        print(gcf,'-dpng',[rootdir 'distribution_test_sources.png']);
    end
    % ---------------------------------------------------------------------
    % Drawing the SVM's predictions of good performance
    for i=1:nalgos
        clf;
        drawBinaryPerformance(out.pilot.Z, out.pythia.Yhat(:,i), strrep(algolabels{i},'_',' '));
        print(gcf,'-dpng',[rootdir 'binary_test_svm_' algolabels{i} '.png']);
    end
    % ---------------------------------------------------------------------
    % Drawing the SVM's recommendations
    clf;
    drawPortfolioSelections(out.pilot.Z, psel, algolabels, 'Predicted best algorithm');
    print(gcf,'-dpng',[rootdir 'distribution_test_svm_portfolio.png']);
end

disp('-------------------------------------------------------------------------');
disp('-> Storing the raw MATLAB results for post-processing and/or debugging.');
save([rootdir 'workspace_test.mat']); % Save the full workspace for debugging
disp(['-> Completed! Elapsed time: ' num2str(toc(startProcess)) 's']);
disp('EOF:SUCCESS');
end