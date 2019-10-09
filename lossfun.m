function out = lossfun(xtr,ytr,xte,yte,param,weight,flag)

% COL1 = True Negatives  (y==0 & yhat==0)
% COL2 = False Negatives (y==1 & yhat==0)
% COL3 = False Positives (y==0 & yhat==1)
% COL4 = True Positives  (y==1 & yhat==1)

out = confusionmat(yte,svmwrap(xtr,ytr,xte,param,weight,flag));
if all(size(out)~=[2 2])
    ctr = sum(bsxfun(@eq,ytr,[1 2]));
    cte = sum(bsxfun(@eq,yte,[1 2]));
    if ctr(1)==0 || cte(1)==0
        out = [0 0; 0 out];
    elseif ctr(2)==0 || cte(2)==0
        out = [out 0; 0 0];
    end
end

end