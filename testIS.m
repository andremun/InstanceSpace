function [out] = testIS(rootdir)
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
printdisclaim('testIS.m');
% -------------------------------------------------------------------------
% Collect all the data from the files
modelfile = [rootdir 'model.mat'];
datafile = [rootdir 'metadata.csv'];
optsfile = [rootdir 'options.json'];
if ~isfile(modelfile) || ~isfile(datafile) || ~isfile(optsfile)
    error('Please place the datafiles in the directory specified.');
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
if model.opts.perf.MaxPerf
    Y(isnan(Y)) = -Inf;
    [bestPerformace,portfolio] = max(Y,[],2);
    if model.opts.perf.AbsPerf
        Ybin = Y>=model.opts.perf.epsilon;
        msg = [msg 'higher than ' num2str(model.opts.perf.epsilon)];
    else
        Ybin = bsxfun(@ge,Y,(1-model.opts.perf.epsilon).*bestPerformace); % One is good, zero is bad
        msg = [msg 'within ' num2str(round(100.*model.opts.perf.epsilon)) '% of the best.'];
    end
else
    Y(isnan(Y)) = Inf;
    [bestPerformace,portfolio] = min(Y,[],2);
    if model.opts.perf.AbsPerf
        Ybin = Y<=model.opts.perf.epsilon;
        msg = [msg 'less than ' num2str(model.opts.perf.epsilon)];
    else
        Ybin = bsxfun(@le,Y,(1+model.opts.perf.epsilon).*bestPerformace);
        msg = [msg 'within ' num2str(round(100.*model.opts.perf.epsilon)) '% of the best.'];
    end
end
disp(msg);
beta = sum(Ybin,2)>model.opts.general.betaThreshold*nalgos;

% ---------------------------------------------------------------------
% Automated pre-processing
if model.opts.bound.flag
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

if model.opts.norm.flag
    % Normalize the data using Box-Cox and Z-transformations
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
Z = X*model.pbldr.A';
% -------------------------------------------------------------------------
% Algorithm selection. Fit a model that would separate the space into
% classes of good and bad performance. 
Yhat = 0.*Ybin;
probs = 0.*Ybin;
for i=1:nalgos
    [Yhat(:,i),~,probs(:,i)] = svmpredict(randi([1 2], ninst,1), Z, model.algosel.svm{i}, '-q');
end
Yhat = Yhat==2; % Make it binary
[mostaccurate,psel] = max(bsxfun(@times,Yhat,model.algosel.precision),[],2);
psel(mostaccurate<=0) = 0;

writeArray2CSV = @(data,colnames,rownames,filename) writetable(array2table(data,'VariableNames',colnames,...
                                                                                'RowNames',rownames),...
                                                               filename,'WriteRowNames',true);
colorscale  = @(data) round(255.*bsxfun(@rdivide, bsxfun(@minus, data, min(data,[],1)), range(data)));
colorscaleg = @(data) round(255.*bsxfun(@rdivide, bsxfun(@minus, data, min(data(:))), range(data(:))));

if opts.outputs.csv
    disp('-------------------------------------------------------------------------');
    disp('-> Writing the data on CSV files for post-processing.');
    % ---------------------------------------------------------------------
    writeArray2CSV(Z, {'z_1','z_2'}, instlabels, [rootdir 'coordinates.csv']);
    writeArray2CSV(Xraw(:, model.featsel.idx), featlabels, instlabels, [rootdir 'feature_raw.csv']);
    writeArray2CSV(X, featlabels, instlabels, [rootdir 'feature_process.csv']);  
    writeArray2CSV(Yraw, algolabels, instlabels, [rootdir 'algorithm_raw.csv']);
    writeArray2CSV(Y, algolabels, instlabels, [rootdir 'algorithm_process.csv']);
    writeArray2CSV(Ybin, algolabels, instlabels, [rootdir 'algorithm_bin.csv']);
    writeArray2CSV(portfolio, {'Best_Algorithm'}, instlabels, [rootdir 'portfolio.csv']);
    writeArray2CSV(Yhat, algolabels, instlabels, [rootdir 'algorithm_svm.csv']);
    writeArray2CSV(psel, {'Best_Algorithm'}, instlabels, [rootdir 'portfolio_svm.csv']);
    if opts.outputs.web
    %   writetable(array2table(parula(256), 'VariableNames', {'R','G','B'}), [rootdir 'color_table.csv']);
        writeArray2CSV(colorscale(Xraw(:,model.featsel.idx)), featlabels, instlabels, [rootdir 'feature_raw_color.csv']);
        writeArray2CSV(colorscale(Yraw), algolabels, instlabels, [rootdir 'algorithm_raw_single_color.csv']);
        writeArray2CSV(colorscale(X), featlabels, instlabels, [rootdir 'feature_process_color.csv']);
        writeArray2CSV(colorscale(Y), algolabels, instlabels, [rootdir 'algorithm_process_single_color.csv']);
        writeArray2CSV(colorscaleg(Yraw), algolabels, instlabels, [rootdir 'algorithm_raw_color.csv']);
        writeArray2CSV(colorscaleg(Y), algolabels, instlabels, [rootdir 'algorithm_process_color.csv']);
    end
end

if opts.outputs.png
    disp('-------------------------------------------------------------------------');
    disp('-> Producing the plots.');
    % ---------------------------------------------------------------------
    for i=1:nfeats
        clf;
        drawScatter((X(:,i)-min(X(:,i)))./range(X(:,i)), Z, strrep(featlabels{i},'_',' '));
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
        drawScatter(Ys(:,i), Z, strrep(algolabels{i},'_',' '));
        print(gcf,'-dpng',[rootdir 'scatter_' algolabels{i} '_absolute.png']);
    end
    % ---------------------------------------------------------------------
    for i=1:nalgos
        clf;
        drawScatter((Y(:,i)-min(Y(:,i)))./range(Y(:,i)), Z, strrep(algolabels{i},'_',' '));
        print(gcf,'-dpng',[rootdir 'scatter_' algolabels{i} '.png']);
    end
    % ---------------------------------------------------------------------
    % Drawing the sources of the instances if available
    if any(issource)
        clf;
        drawSources(Z, S);
        print(gcf,'-dpng',[rootdir 'sources.png']);
    end
    % ---------------------------------------------------------------------
    % Drawing the SVM's predictions of good performance
    for i=1:nalgos
        clf;
        drawSVMPredictions(Z, Yhat(:,i), strrep(algolabels{i},'_',' '));
        print(gcf,'-dpng',[rootdir 'svm_' algolabels{i} '.png']);
    end
    % ---------------------------------------------------------------------
    % Drawing the SVM's recommendations
    clf;
    drawSVMPortfolioSelections(Z, psel, algolabels);
    print(gcf,'-dpng',[rootdir 'svm_portfolio.png']);
end

disp('-------------------------------------------------------------------------');
disp('-> Storing the raw MATLAB results for post-processing and/or debugging.');
save([rootdir 'workspace_test.mat']); % Save the full workspace for debugging
disp(['-> Completed! Elapsed time: ' num2str(toc(startProcess)) 's']);
disp('EOF:SUCCESS');

end