function [X, out] = checkCorrelation(X,Y,opts)

[~,nfeats] = size(X);
if nfeats>2
    disp('-> Checking for feature correlation with performance.');
    [out.rho,out.p] = corr(X,Y,'rows','pairwise');
    rho = out.rho;
    rho(isnan(rho) | (out.p>0.05)) = 0;
    [~,row] = sort(abs(rho),1,'descend');
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
    disp('-> There are less than 3 features to do selection (correlation). Skipping.')
    out.p = ones(nfeats,size(Y,2));
    out.rho = NaN.*out.p;
    out.selvars = true(1,nfeats);
end

end