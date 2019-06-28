function [X, out] = checkCorrelation(X,Y,opts)

[~,nfeats] = size(X);
if opts.flag && nfeats>2
    disp('-> Checking for feature correlation with performance.');
    out.rho = corr(X,Y,'rows','pairwise');
    out.rho(isnan(out.rho)) = 0;
    [~,row] = sort(abs(out.rho),1,'descend');
    out.selvars = false(1,nfeats);
    testTreshold = false;
    while sum(out.selvars)<3
        if testTreshold
            warning('Feature selection using correlation was too strict. The threshold value was increased automatically by one.')
        end
        out.selvars(unique(row(1:min(opts.threshold,nfeats),:))) = true;
        opts.threshold = opts.threshold + 1;
        testTreshold = true;
    end
    X = X(:,out.selvars);
    disp(['-> Keeping ' num2str(size(X,2)) ' out of ' num2str(nfeats) ' features (correlation).']);
else
    out.rho = NaN.*ones(nfeats,size(Y,2));
    out.selvars = true(1,nfeats);
end

end