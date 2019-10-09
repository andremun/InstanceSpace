function out = fitoracle(Z,Ybin,weight,opts)

% global params
Ybin = double(Ybin)+1;
nalgos = size(Ybin,2);
out.paramgrid = sortrows(2.^(opts.maxcvgrid.*lhsdesign(opts.cvgrid,2) + ...
                             opts.mincvgrid));  % Cross-validation grid
out.cvcmat = NaN.*ones(opts.cvgrid,4,nalgos);
out.paramidx = NaN.*ones(1,nalgos);
out.precision = NaN.*ones(1,nalgos);
out.recall = NaN.*ones(1,nalgos);
out.accuracy = NaN.*ones(1,nalgos);

for i=1:nalgos
    cp = cvpartition(Ybin(:,i),'Kfold',opts.cvfolds);
    for j=1:opts.cvgrid
        f = @(xtr,ytr,xte,yte) lossfun(xtr,ytr,xte,yte,out.paramgrid(j,:),weight(i),false);
        out.cvcmat(j,:,i) = sum(crossval(f, Z, Ybin(:,i),'partition', cp));
    end
    tn = out.cvcmat(:,1,i);
    fp = out.cvcmat(:,3,i);
    fn = out.cvcmat(:,2,i);
    tp = out.cvcmat(:,4,i);
    [out.precision(i),best] = max(tp./(tp+fp));
    out.paramidx(i) = best;
    out.recall(i) = tp(best)./(tp(best)+fn(best));
    out.accuracy(i) = (tp(best)+tn(best))./sum(out.cvcmat(best,:,i));
    disp(['    -> ' num2str(i) ' out of ' num2str(nalgos) ' models have been fitted.']);
end
disp(['    -> Completed - Average cross validation precision is: ' ...
      num2str(round(100.*mean(out.precision),1)) '%']);

out.Yhat = 0.*Ybin;
out.probs = 0.*Ybin;
out.svm = cell(1,nalgos);
out.svmparams = zeros(nalgos,2);
for i=1:nalgos
    out.svmparams(i,:) = out.paramgrid(out.paramidx(i),:);
    [aux, out.svm{i}] = svmwrap(Z, Ybin(:,i), Z, out.svmparams(i,:), weight(i), true);
    out.Yhat(:,i)  = aux(:,1);
    out.probs(:,i) = aux(:,2);
end
out.Yhat = out.Yhat==2; % Make it binary
% We assume that the most precise SVM (as per CV-Error) is the most
% reliable.
% [mostaccurate,out.psel] = max(bsxfun(@times,out.Yhat,1-out.modelerr),[],2);
[mostprecise,out.psel] = max(bsxfun(@times,out.Yhat,out.precision),[],2);
out.pselfull = out.psel;
out.psel(mostprecise<=0) = 0;
[~,betterdefault] = max(mean(Ybin));
out.pselfull(mostprecise<=0) = betterdefault;

end