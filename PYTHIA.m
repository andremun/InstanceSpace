function out = PYTHIA(Z, Y, Ybin, Ybest, algolabels, opts)
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
[Znorm,out.mu,out.sigma] = zscore(Z);
[ninst,nalgos] = size(Ybin);
out.cp = cell(1,nalgos);
out.svm = cell(1,nalgos);
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
params = NaN.*ones(nalgos,2);
if opts.ispolykrnl
    KernelFcn = 'polynomial';
else
    if ninst>1e3
        disp('  -> For datasets larger than 1K Instances, PYTHIA works better with a Polynomial kernel.');
        disp('  -> Consider changing the kernel if the results are unsatisfactory.');
        disp('-------------------------------------------------------------------------');
    end
    KernelFcn = 'gaussian';
end
disp(['  -> PYTHIA is using a ' KernelFcn ' kernel ']);
disp('-------------------------------------------------------------------------');
if opts.uselibsvm
    disp('  -> Using LIBSVM''s libraries.');
    if precalcparams
        disp('  -> Using pre-calculated hyper-parameters for the SVM.');
        params = opts.params;
    else
        disp('  -> Search on a latin hyper-cube design will be used for parameter hyper-tunning.');
    end
else
    disp('  -> Using MATLAB''s SVM libraries.');
    if precalcparams
        disp('  -> Using pre-calculated hyper-parameters for the SVM.');
        params = opts.params;
    else
        disp('  -> Bayesian Optimization will be used for parameter hyper-tunning.');
    end
    disp('-------------------------------------------------------------------------');
    if opts.useweights
        disp('  -> PYTHIA is using cost-sensitive classification');
        out.W = abs(Y-nanmean(Y(:)));
        out.W(out.W==0) = min(out.W(out.W~=0));
        out.W(isnan(out.W)) = max(out.W(~isnan(out.W)));
        Waux = out.W;
    else
        disp('  -> PYTHIA is not using cost-sensitive classification');
        Waux = ones(ninst,nalgos);
    end
end
disp('-------------------------------------------------------------------------');
disp(['  -> Using a ' num2str(opts.cvfolds) ...
      '-fold stratified cross-validation experiment to evaluate the SVMs.']);
disp('-------------------------------------------------------------------------');
disp('  -> Training has started. PYTHIA may take a while to complete...');
t = tic;
for i=1:nalgos
    tic;
    state = rng;
    rng('default');
    out.cp{i} = cvpartition(Ybin(:,i),'Kfold',opts.cvfolds,'Stratify',true);
    if opts.uselibsvm
        [out.svm{i},out.Ysub(:,i),out.Pr0sub(:,i),out.Yhat(:,i),...
            out.Pr0hat(:,i),out.boxcosnt(i),out.kscale(i)] = fitlibsvm(Znorm,Ybin(:,i),...
                                                                       out.cp{i},KernelFcn,...
                                                                       params(i,:));
    else
        [out.svm{i},out.Ysub(:,i),out.Pr0sub(:,i),out.Yhat(:,i),...
            out.Pr0hat(:,i),out.boxcosnt(i),out.kscale(i)] = fitmatsvm(Znorm,Ybin(:,i),...
                                                                       Waux(:,i),out.cp{i},...
                                                                       KernelFcn,params(i,:));
    end
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
                        'BoxConstraint';
                        'KernelScale'};
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
function [svm,Ysub,Psub,Yhat,Phat,C,g] = fitlibsvm(Z,Ybin,cp,k,params)

ninst = size(Z,1);
maxgrid =   4;
mingrid = -10;
if any(isnan(params))
    rng('default');
    nvals = 30;
    paramgrid = sortrows(2.^((maxgrid-mingrid).*lhsdesign(nvals,2) + mingrid));
else
    nvals = 1;
    paramgrid = params;
end
Ybin = double(Ybin)+1;
Ysub = zeros(ninst,nvals);
Psub = zeros(ninst,nvals);

if strcmp(k,'polynomial')
    k = 1;
else
    k = 2;
end

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

for jj=1:cp.NumTestSets
    idx = cp.training(jj);
    Ztrain = Z(idx,:);
    Ytrain = Ybin(idx);
    Ztest = Z(~idx,:);
    Ytest = Ybin(~idx);
    Yaux = zeros(sum(~idx),nvals);
    Paux = zeros(sum(~idx),nvals);
    parfor (ii=1:nvals,nworkers)
        cparams = paramgrid(ii,:);
        prior = mean(bsxfun(@eq,Ytrain,[1 2]));
        command = ['-s 0 -t ' num2str(k) ' -q -b 1 -c ' num2str(cparams(1)) ...
                   ' -g ' num2str(cparams(2)) ' -w1 1 -w2 ' num2str(prior(1)./prior(2),4)];
        rng('default');
        svm = svmtrain(Ytrain, Ztrain, command); %#ok<SVMTRAIN>
        [Yaux(:,ii),~,Paux(:,ii)] = svmpredict(Ytest, Ztest, svm, '-q');
    end
    for ii=1:nvals
        Ysub(~idx,ii) = Yaux(:,ii);
        Psub(~idx,ii) = Paux(:,ii);
    end
end
[~,idx] = min(mean(bsxfun(@ne,Ysub,Ybin),1));
Ysub = Ysub(:,idx)==2;
Psub = Psub(:,idx);

C = paramgrid(idx,1);
g = paramgrid(idx,2);
prior = mean(bsxfun(@eq,Ybin,[1 2]));
command = ['-s 0 -t ' num2str(k) ' -q -b 1 -c ' num2str(C) ' -g ' num2str(g) ...
           ' -w1 1 -w2 ' num2str(prior(1)./prior(2),4)];
rng('default');
svm = svmtrain(Ybin, Z, command); %#ok<SVMTRAIN>
[Yhat,~,Phat] = svmpredict(Ybin, Z, svm, '-q');
Yhat = Yhat==2;

end
% =========================================================================
function [svm,Ysub,Psub,Yhat,Phat,C,g] = fitmatsvm(Z,Ybin,W,cp,k,params)

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

if any(isnan(params))
    rng('default');
    hypparams = hyperparameters('fitcsvm',Z,Ybin);
    hypparams = hypparams(1:2);
    hypparams(1).Range = 2.^[-10,4];
    hypparams(2).Range = hypparams(1).Range;
    svm = fitcsvm(Z,Ybin,'Standardize',false,...
                         'Weights',W,...
                         'CacheSize','maximal',...
                         'RemoveDuplicates',true,...
                         'KernelFunction',k,...
                         'OptimizeHyperparameters',hypparams,...
                         'HyperparameterOptimizationOptions',...
                         struct('CVPartition',cp,...
                                'Verbose',0,...
                                'AcquisitionFunctionName','probability-of-improvement',...
                                'ShowPlots',false,...
                                'UseParallel',nworkers~=0));
    svm = fitSVMPosterior(svm);
    C = svm.HyperparameterOptimizationResults.bestPoint{1,1};
    g = svm.HyperparameterOptimizationResults.bestPoint{1,2};
    [Ysub,aux] = svm.resubPredict;
    Psub = aux(:,1);
    [Yhat,aux] = svm.predict(Z);
    Phat = aux(:,1);
else
    C = params(1);
    g = params(2);
    rng('default');
    svm = fitcsvm(Z,Ybin,'Standardize',false,...
                         'Weights',W,...
                         'CacheSize','maximal',...
                         'RemoveDuplicates',true,...
                         'KernelFunction',k,...
                         'CVPartition',cp,...
                         'BoxConstraint',C,...
                         'KernelScale',g);
    svm = fitSVMPosterior(svm);
    [Ysub,aux] = svm.kfoldPredict;
    Psub = aux(:,1);
    rng('default');
    svm = fitcsvm(Z,Ybin,'Standardize',false,...
                         'Weights',W,...
                         'CacheSize','maximal',...
                         'RemoveDuplicates',true,...
                         'KernelFunction',k,...
                         'BoxConstraint',C,...
                         'KernelScale',g);
	svm = fitSVMPosterior(svm);
    [Yhat,aux] = svm.predict(Z);
    Phat = aux(:,1);
end

end
% =========================================================================
