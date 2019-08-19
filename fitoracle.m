function out = fitoracle(Z,Ybin,opts)

% global params
Ybin = double(Ybin)+1;
nalgos = size(Ybin,2);
out.paramgrid = sortrows(2.^(opts.maxcvgrid.*lhsdesign(opts.cvgrid,2) + ...
                             opts.mincvgrid));  % Cross-validation grid
out.cvmcr = NaN.*ones(opts.cvgrid,nalgos);
out.paramidx = NaN.*ones(1,nalgos);
out.modelerr = NaN.*ones(1,nalgos);

for i=1:nalgos
    for j=1:opts.cvgrid
        out.cvmcr(j,i) = crossval('mcr', Z, Ybin(:,i),...
                                  'Kfold', opts.cvfolds,...
                                  'Options', statset('UseParallel',true),...
                                  'Predfun',@(xtrain,ytrain,xtest) svmwrap(xtrain,...
                                                                           ytrain,...
                                                                           xtest,...
                                                                           out.paramgrid(j,:),...
                                                                           false));
    end
    [out.modelerr(i),out.paramidx(i)] = min(out.cvmcr(:,i));
    disp(['    ->  ' num2str(i) ' out of ' num2str(nalgos) ' models have been fitted.']);
end

out.Yhat = 0.*Ybin;
out.probs = 0.*Ybin;
out.svm = cell(1,nalgos);
for i=1:nalgos
    [aux, out.svm{i}] = svmwrap(Z, Ybin(:,i), Z, out.paramgrid(out.paramidx(i),:), true);
    out.Yhat(:,i)  = aux(:,1);
    out.probs(:,i) = aux(:,2);
end
out.Yhat = out.Yhat==2; % Make it binary
% [~,out.psel] = max(out.probs,[],2); % Determine which one to suggest
% - Different strategy to make the SVM recommendation. Use the SVM which
% has the lowest error.
[mostaccurate,out.psel] = max(bsxfun(@times,out.Yhat,1-out.modelerr),[],2);
out.psel(mostaccurate<=0) = 0;