function out = PYTHIA2(Z, Y, Ybin, Ybest, algolabels, opts)
% -------------------------------------------------------------------------
% PYTHIA.m
% -------------------------------------------------------------------------
%
% By: Mario Andres Munoz Acosta
%     School of Mathematics and Statistics
%     The University of Melbourne
%     Australia
%     2020
%
% -------------------------------------------------------------------------

disp('  -> Initializing PYTHIA.');
[ninst,nalgos] = size(Ybin);
out.cp = cell(1,nalgos);
out.knn = cell(1,nalgos);
out.cvcmat = zeros(nalgos,4);
out.Ysub = false & Ybin;
out.Yhat = false & Ybin;
out.Pr0sub = 0.*Ybin;
out.Pr0hat = 0.*Ybin;
out.boxcosnt = zeros(1,nalgos);
out.kscale = out.boxcosnt;
disp('-------------------------------------------------------------------------');
precalcparams = isfield(opts,'params') && isnumeric(opts.params) && ...
                size(opts.params,1)==nalgos && size(opts.params,2)==2;
disp('  -> Using MATLAB''s KNN libraries.');
if precalcparams
    disp('  -> Using pre-calculated hyper-parameters for the SVM.');
    params = opts.params;
else
    disp('  -> Bayesian Optimization will be used for parameter hyper-tunning.');
    params = NaN.*ones(nalgos,2);
end
disp('-------------------------------------------------------------------------');
if opts.useweights
    disp('  -> PYTHIA is using weights for cost-sensitive classification [experimental]');
    out.W = abs(Y-nanmean(Y(:)));
    out.W(out.W==0) = min(out.W(out.W~=0));
    out.W(isnan(out.W)) = max(out.W(~isnan(out.W)));
    Waux = out.W;
else
    disp('  -> PYTHIA is not using cost-sensitive classification');
    Waux = ones(ninst,nalgos);
end

disp('-------------------------------------------------------------------------');
disp(['  -> Using a ' num2str(opts.cvfolds) ...
      '-fold stratified cross-validation experiment to evaluate the models.']);
disp('-------------------------------------------------------------------------');
disp('  -> Training has started. PYTHIA may take a while to complete...');
t = tic;
for i=1:nalgos
    tic;
    state = rng;
    rng('default');
    out.cp{i} = cvpartition(Ybin(:,i),'Kfold',opts.cvfolds,'Stratify',true);
    [out.knn{i},out.Ysub(:,i),out.Pr0sub(:,i),out.Yhat(:,i),...
        out.Pr0hat(:,i),out.boxcosnt(i),out.kscale(i)] = fitmatknn(Z,Ybin(:,i),...
                                                                   Waux(:,i),out.cp{i},...
                                                                   params(i,:));
    rng(state);
    aux = confusionmat(Ybin(:,i),out.Ysub(:,i));
    if numel(aux)~=4
        caux = aux;
        aux = zeros(2);
        if all(Ybin(:,i)==0)
            if all(out.Ysub(:,i)==0)
                aux(1,1) = caux;
            elseif all(out.Ysub(:,i)==1)
                aux(2,1) = caux;
            end
        elseif all(Ybin(:,i)==1)
            if all(out.Ysub(:,i)==0)
                aux(1,2) = caux;
            elseif all(out.Ysub(:,i)==1)
                aux(2,2) = caux;
            end
        end
    end
	out.cvcmat(i,:) = aux(:);
    if i==nalgos
        disp(['    -> PYTHIA has trained a model for ''' algolabels{i}, ...
              ''', there are no models left to train.']);
    elseif i==nalgos-1
        disp(['    -> PYTHIA has trained a model for ''' algolabels{i}, ...
              ''', there is 1 model left to train.']);
    else
        disp(['    -> PYTHIA has trained a model for ''' algolabels{i}, ...
              ''', there are ' num2str(nalgos-i) ' models left to train.']);
    end
    disp(['      -> Elapsed time: ' num2str(toc,'%.2f\n') 's']);
end
tn = out.cvcmat(:,1);
fp = out.cvcmat(:,3);
fn = out.cvcmat(:,2);
tp = out.cvcmat(:,4);
out.precision = tp./(tp+fp);
out.recall = tp./(tp+fn);
out.accuracy = (tp+tn)./ninst;
disp('-------------------------------------------------------------------------');
disp('  -> PYTHIA has completed training the models.');
disp(['  -> The average cross validated precision is: ' ...
      num2str(round(100.*mean(out.precision),1)) '%']);
disp(['  -> The average cross validated accuracy is: ' ...
      num2str(round(100.*mean(out.accuracy),1)) '%']);
  disp(['      -> Elapsed time: ' num2str(toc(t),'%.2f\n') 's']);
disp('-------------------------------------------------------------------------');
% We assume that the most precise SVM (as per CV-Precision) is the most
% reliable.
if nalgos>1
    [best,out.selection0] = max(bsxfun(@times,out.Yhat,out.precision'),[],2);
else
    best = out.Yhat;
    out.selection0 = out.Yhat;
end
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
out.summary = cell(nalgos+3, 11);
out.summary{1,1} = 'Algorithms ';
out.summary(2:end-2, 1) = algolabels;
out.summary(end-1:end, 1) = {'Oracle','Selector'};
out.summary(1, 2:11) = {'Avg_Perf_all_instances';
                        'Std_Perf_all_instances';
                        'Probability_of_good';
                        'Avg_Perf_selected_instances';
                        'Std_Perf_selected_instances';
                        'CV_model_accuracy';
                        'CV_model_precision';
                        'CV_model_recall';
                        'NumNeighbours';
                        'Distance'};
out.summary(2:end, 2) = num2cell(round([avgperf nanmean(Ybest) nanmean(Yfull(:))],3));
out.summary(2:end, 3) = num2cell(round([stdperf nanstd(Ybest) nanstd(Yfull(:))],3));
out.summary(2:end, 4) = num2cell(round([mean(Ybin) 1 pgood],3));
out.summary(2:end, 5) = num2cell(round([nanmean(Ysvms) NaN nanmean(Y(:))],3));
out.summary(2:end, 6) = num2cell(round([nanstd(Ysvms) NaN nanstd(Y(:))],3));
out.summary(2:end, 7) = num2cell(round(100.*[out.accuracy' NaN NaN],1));
out.summary(2:end, 8) = num2cell(round(100.*[out.precision' NaN precisionsel],1));
out.summary(2:end, 9) = num2cell(round(100.*[out.recall' NaN recallsel],1));
out.summary(2:end-2, 10) = num2cell(round(out.boxcosnt,3));
out.summary(2:end-2, 11) = num2cell(round(out.kscale,3));
out.summary(cellfun(@(x) all(isnan(x)),out.summary)) = {[]}; % Clean up. Not really needed
disp('  -> PYTHIA has completed! Performance of the models:');
disp(' ');
disp(out.summary);

end
% =========================================================================
% SUBFUNCTIONS
% =========================================================================
function [knn,Ysub,Psub,Yhat,Phat,NumNeighbors,Distance] = fitmatknn(Z,Ybin,W,cp,params)

if exist('gcp','file')==2
    mypool = gcp('nocreate');
    if ~isempty(mypool)
        nworkers = mypool.NumWorkers;
    else
        nworkers = 0;
    end
else
    nworkers = 0;
end

hypparams = hyperparameters('fitcknn',Z,Ybin);

if any(isnan(params))
    rng('default');
    knn = fitcknn(Z,Ybin,'Standardize',false,...
                         'Weights',W,...
                         'OptimizeHyperparameters',hypparams,...
                         'HyperparameterOptimizationOptions',...
                         struct('CVPartition',cp,...
                                'Verbose',0,...
                                'AcquisitionFunctionName','probability-of-improvement',...
                                'ShowPlots',false,...
                                'UseParallel',nworkers~=0));
    NumNeighbors = knn.HyperparameterOptimizationResults.bestPoint{1,1};
    Distance = knn.HyperparameterOptimizationResults.bestPoint{1,2};
    [Ysub,aux] = knn.resubPredict;
    Psub = aux(:,1);
    [Yhat,aux] = knn.predict(Z);
    Phat = aux(:,1);
else
    NumNeighbors = round(params(1));
    Distance = hypparams(2).Range(round(params(2)));
    rng('default');
    knn = fitcknn(Z,Ybin,'Standardize',false,...
                         'Weights',W,...
                         'CVPartition',cp,...
                         'NumNeighbors',NumNeighbors,...
                         'Distance',Distance);
    [Ysub,aux] = knn.kfoldPredict;
    Psub = aux(:,1);
    [Yhat,aux] = knn.predict(Z);
    Phat = aux(:,1);
end

end
% =========================================================================