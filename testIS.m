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
datafile = [rootdir 'metadata_test.csv'];
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
out.data.instlabels = Xbar{:,isname};
if isnumeric(out.data.instlabels)
    out.data.instlabels = num2cell(out.data.instlabels);
    out.data.instlabels = cellfun(@(x) num2str(x),out.data.instlabels,'UniformOutput',false);
end
if any(issource)
    out.data.S = categorical(Xbar{:,issource});
end
out.data.X = Xbar{:,isfeat};
out.data.Y = Xbar{:,isalgo};
nalgos = size(out.data.Y,2);
% -------------------------------------------------------------------------
% Storing the raw data for further processing, e.g., graphs
out.data.Xraw = out.data.X;
out.data.Yraw = out.data.Y;
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
    Yaux = out.data.Y;
    Yaux(isnan(Yaux)) = -Inf;
    [rankPerf,rankAlgo] = sort(Yaux,2,'descend');
    out.data.bestPerformace = rankPerf(:,1);
    out.data.P = rankAlgo(:,1);
    if model.opts.perf.AbsPerf
        out.data.Ybin = out.data.Y>=model.opts.perf.epsilon;
        msg = [msg 'higher than ' num2str(model.opts.perf.epsilon)];
    else
        out.data.Ybin = bsxfun(@ge,out.data.Y,(1-model.opts.perf.epsilon).*out.data.bestPerformace); % One is good, zero is bad
        msg = [msg 'within ' num2str(round(100.*model.opts.perf.epsilon)) '% of the best.'];
    end
else
    Yaux = out.data.Y;
    Yaux(isnan(Yaux)) = Inf;
    [rankPerf,rankAlgo] = sort(Yaux,2,'ascend');
    out.data.bestPerformace = rankPerf(:,1);
    out.data.P = rankAlgo(:,1);
    if model.opts.perf.AbsPerf
        out.data.Ybin = out.data.Y<=model.opts.perf.epsilon;
        msg = [msg 'less than ' num2str(model.opts.perf.epsilon)];
    else
        out.data.Ybin = bsxfun(@le,out.data.Y,(1+model.opts.perf.epsilon).*out.data.bestPerformace);
        msg = [msg 'within ' num2str(round(100.*model.opts.perf.epsilon)) '% of the best.'];
    end
end
disp(msg);
out.data.numGoodAlgos = sum(out.data.Ybin,2);
out.data.beta = out.data.numGoodAlgos>model.opts.general.betaThreshold*nalgos;

% ---------------------------------------------------------------------
% Automated pre-processing
if model.opts.auto.preproc && model.opts.bound.flag
    disp('-------------------------------------------------------------------------');
    disp('-> Auto-pre-processing. Bounding outliers, scaling and normalizing the data.');
    % Eliminate extreme outliers, i.e., any point that exceedes 5 times the
    % inter quantile range, by bounding them to that value.
    disp('-> Removing extreme outliers from the feature values.');
    himask = bsxfun(@gt,out.data.X,model.bound.hibound);
    lomask = bsxfun(@lt,out.data.X,model.bound.lobound);
    out.data.X = out.data.X.*~(himask | lomask) + bsxfun(@times,himask,model.bound.hibound) + ...
                                                  bsxfun(@times,lomask,model.bound.lobound);
end

if model.opts.auto.preproc && model.opts.norm.flag
    % Normalize the data using Box-Cox and out.pilot.Z-transformations
    disp('-> Auto-normalizing the data.');
    out.data.X = bsxfun(@minus,out.data.X,model.norm.minX)+1;
    out.data.X = bsxfun(@rdivide,bsxfun(@power,out.data.X,model.norm.lambdaX)-1,model.norm.lambdaX);
    out.data.X = bsxfun(@rdivide,bsxfun(@minus,out.data.X,model.norm.muX),model.norm.sigmaX);
    
    out.data.Y(out.data.Y==0) = eps; % Assumes that out.data.Y is always positive and higher than 1e-16
    out.data.Y = bsxfun(@rdivide,bsxfun(@power,out.data.Y,model.norm.lambdaY)-1,model.norm.lambdaY);
    out.data.Y = bsxfun(@rdivide,bsxfun(@minus,out.data.Y,model.norm.muY),model.norm.sigmaY);
end
% ---------------------------------------------------------------------
% This is the final subset of features.
out.featsel.idx = model.featsel.idx;
out.data.X = out.data.X(:,out.featsel.idx);
out.data.featlabels = strrep(varlabels(isfeat),'feature_','');
out.data.featlabels = out.data.featlabels(model.featsel.idx);
out.data.algolabels = strrep(varlabels(isalgo),'algo_','');
% ---------------------------------------------------------------------
%  Calculate the two dimensional projection using the PBLDR algorithm
%  (Munoz et al. Mach Learn 2018)
out.pilot.Z = out.data.X*model.pilot.A';
% -------------------------------------------------------------------------
% Algorithm selection. Fit a model that would separate the space into
% classes of good and bad performance. 
out.pythia = PYTHIAtest(model.pythia, out.pilot.Z, out.data.Yraw, ...
                        out.data.Ybin, out.data.bestPerformace, ...
                        out.data.algolabels);
% -------------------------------------------------------------------------
% Validating the footprints
if opts.trace.usesim
    out.trace = TRACEtest(model.trace, out.pilot.Z, out.data.Ybin, ...
                          out.pythia.selection0, out.data.beta, ...
                          out.data.algolabels);
else
    out.trace = TRACEtest(model.trace, out.pilot.Z, out.data.Ybin, ...
                          out.data.P, out.data.beta, ...
                          out.data.algolabels);
end

out.opts = opts;
% -------------------------------------------------------------------------
% Writing the results
if opts.outputs.csv
    scriptcsv(out,rootdir);
    if opts.outputs.web
        scriptweb(out,rootdir);
    end
end

if opts.outputs.png
    scriptpng(out,rootdir);
end

disp('-------------------------------------------------------------------------');
disp('-> Storing the raw MATLAB results for post-processing and/or debugging.');
save([rootdir 'workspace_test.mat']); % Save the full workspace for debugging
disp(['-> Completed! Elapsed time: ' num2str(toc(startProcess)) 's']);
disp('EOF:SUCCESS');
end
% =========================================================================
% SUBFUNCTIONS
% =========================================================================
function out = PYTHIAtest(model, Z, Y, Ybin, Ybest, algolabels)

nalgos = size(Ybin,2);
out.Yhat = false(size(Ybin));
out.Pr0hat = 0.*Ybin;
out.cvcmat = zeros(nalgos,4);
for i=1:nalgos
    [out.Yhat(:,i),aux] = model.svm{i}.predict(Z);
    out.Pr0hat(:,i) = aux(:,1);
    aux = confusionmat(Ybin(:,i),out.Yhat(:,i));
    out.cvcmat(i,:) = aux(:);
end
tn = out.cvcmat(:,1);
fp = out.cvcmat(:,3);
fn = out.cvcmat(:,2);
tp = out.cvcmat(:,4);
out.precision = tp./(tp+fp);
out.recall = tp./(tp+fn);
out.accuracy = (tp+tn)./sum(out.cvcmat(1,:));

[best,out.selection0] = max(bsxfun(@times,out.Yhat,model.precision'),[],2);
[~,default] = max(mean(Ybin));
out.selection1 = out.selection0;
out.selection0(best<=0) = 0;
out.selection1(best<=0) = default;

sel0 = bsxfun(@eq,out.selection0,1:nalgos);
sel1 = bsxfun(@eq,out.selection1,1:nalgos);
avgperf = nanmean(Y);
stdperf = nanstd(Y);
Yfull = Y;
Ysvms = Y;
Y(~sel0) = NaN;
Yfull(~sel1) = NaN;
Ysvms(~out.Yhat) = NaN;

pgood = mean(any( Ybin & sel1,2));
fb = sum(any( Ybin & ~sel0,2));
fg = sum(any(~Ybin &  sel0,2));
tg = sum(any( Ybin &  sel0,2));
precisionsel = tg./(tg+fg);
recallsel = tg./(tg+fb);

disp('  -> PYTHIA is preparing the summary table.');
out.summary = cell(nalgos+3, 9);
out.summary{1,1} = 'Algorithms ';
out.summary(2:end-2, 1) = algolabels;
out.summary(end-1:end, 1) = {'Oracle','Selector'};
out.summary(1, 2:9) = {'Avg_Perf_all_instances';
                       'Std_Perf_all_instances';
                       'Probability_of_good';
                       'Avg_Perf_selected_instances';
                       'Std_Perf_selected_instances';
                       'CV_model_accuracy';
                       'CV_model_precision';
                       'CV_model_recall'};
out.summary(2:end, 2) = num2cell(round([avgperf nanmean(Ybest) nanmean(Yfull(:))],3));
out.summary(2:end, 3) = num2cell(round([stdperf nanstd(Ybest) nanstd(Yfull(:))],3));
out.summary(2:end, 4) = num2cell(round([mean(Ybin) 1 pgood],3));
out.summary(2:end, 5) = num2cell(round([nanmean(Ysvms) NaN nanmean(Y(:))],3));
out.summary(2:end, 6) = num2cell(round([nanstd(Ysvms) NaN nanstd(Y(:))],3));
out.summary(2:end, 7) = num2cell(round(100.*[out.accuracy' NaN NaN],1));
out.summary(2:end, 8) = num2cell(round(100.*[out.precision' NaN precisionsel],1));
out.summary(2:end, 9) = num2cell(round(100.*[out.recall' NaN recallsel],1));
out.summary(cellfun(@(x) all(isnan(x)),out.summary)) = {[]}; % Clean up. Not really needed
disp('  -> PYTHIA has completed! Performance of the models:');
disp(' ');
disp(out.summary);

end
% =========================================================================
function model = TRACEtest(model, Z, Ybin, P, beta, algolabels)

nalgos = size(Ybin,2);
disp('-------------------------------------------------------------------------');
disp('  -> TRACE is calculating the algorithm footprints.');
model.test.best = zeros(nalgos,5);
model.test.good = zeros(nalgos,5);
model.test.bad = zeros(nalgos,5);
% Use the actual data to calculate the footprints
for i=1:nalgos
    model.test.best(i,:) = TRACEtestsummary(model.best{i}, Z,  P==i, model.space.area, model.space.density);
    model.test.good(i,:) = TRACEtestsummary(model.good{i}, Z,  Ybin(:,i), model.space.area, model.space.density);
    model.test.bad(i,:)  = TRACEtestsummary(model.bad{i},  Z, ~Ybin(:,i), model.space.area, model.space.density);
end

% -------------------------------------------------------------------------
% Beta hard footprints. First step is to calculate them.
disp('-------------------------------------------------------------------------');
disp('  -> TRACE is calculating the beta-footprints.');
model.test.easy = TRACEtestsummary(model.easy, Z,  beta, model.space.area, model.space.density);
model.test.hard = TRACEtestsummary(model.hard, Z, ~beta, model.space.area, model.space.density);
% -------------------------------------------------------------------------
% Calculating performance
disp('-------------------------------------------------------------------------');
disp('  -> TRACE is preparing the summary table.');
model.summary = cell(nalgos+1,11);
model.summary(1,2:end) = {'Area_Good',...
                          'Area_Good_Normalized',...
                          'Density_Good',...
                          'Density_Good_Normalized',...
                          'Purity_Good',...
                          'Area_Best',...
                          'Area_Best_Normalized',...
                          'Density_Best',...
                          'Density_Best_Normalized',...
                          'Purity_Best'};
model.summary(2:end,1) = algolabels;
model.summary(2:end,2:end) = num2cell([model.test.good model.test.best]);

disp('  -> TRACE has completed. Footprint analysis results:');
disp(' ');
disp(model.summary);

end
% =========================================================================
function out = TRACEtestsummary(footprint, Z, Ybin, spaceArea, spaceDensity)
% 
if isempty(footprint.polygon)
    out = zeros(5,1);
else
    elements = sum(isinterior(footprint.polygon, Z));
    goodElements = sum(isinterior(footprint.polygon, Z(Ybin,:)));
    density = elements./footprint.area;
    purity = goodElements./elements;
    
    out = [footprint.area,...
           footprint.area/spaceArea,...
           density,...
           density/spaceDensity,...
           purity];
end
out(isnan(out)) = 0;
end
% =========================================================================