function [X, out] = checkCorrelation(X,Y,opts)

[~,nfeats] = size(X);
if opts.flag
    disp('-> Checking for feature correlation with performance.');
    out.rho = corr(X,Y,'rows','pairwise');
    [~,row] = sort(abs(out.rho),1,'descend');
    out.selvars = false(1,nfeats);
    out.selvars(unique(row(1:min(opts.threshold,nfeats),:))) = true;
    X = X(:,out.selvars);
    disp(['-> Keeping ' num2str(size(X,2)) ' out of ' num2str(nfeats) ' features (correlation).']);
else
    out.rho = NaN.*ones(nfeats,size(Y,2));
    out.selvars = true(1,nfeats);
end

end