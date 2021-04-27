function out = PYTHIAtest(model, Z, Y, Ybin, Ybest, algolabels)

Z = (Z-model.mu)./model.sigma;
nalgos = length(model.svm);
Y = Y(:,1:nalgos);
Ybin = Ybin(:,1:nalgos);
out.Yhat = false(size(Ybin));
out.Pr0hat = 0.*Ybin;
out.cvcmat = zeros(nalgos,4);
for ii=1:nalgos
    if isstruct(model.svm{ii})
        Yin = double(Ybin(:,ii))+1;
        [aux,~,out.Pr0hat(:,ii)] = svmpredict(Yin, Z, model.svm{ii}, '-q');
        out.Yhat(:,ii) = aux==2;
    elseif isa(model.svm{ii},'ClassificationSVM')
        [out.Yhat(:,ii),aux] = model.svm{ii}.predict(Z);
        out.Pr0hat(:,ii) = aux(:,1);
    else
        disp('There is no model for this algorithm');
        out.Yhat(:,ii) = NaN;
        out.Pr0hat(:,ii) = NaN;
    end
    aux = confusionmat(Ybin(:,ii),out.Yhat(:,ii));
    out.cvcmat(ii,:) = aux(:);
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
out.summary(2:end-2, 1) = algolabels(1:nalgos);
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