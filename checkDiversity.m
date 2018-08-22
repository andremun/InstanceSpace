function [X, out] = checkDiversity(X, opts)

[ninst,nfeats] = size(X);
if opts.flag
    disp('-> Checking for feature diversity.');
    out.pctage = zeros(1,nfeats);
    for i=1:nfeats
        out.pctage(i) = length(unique(X(:,i)))./ninst;
    end
    out.selvars = out.pctage>opts.threshold;
    X = X(:,out.selvars);
    disp(['-> Keeping ' num2str(size(X,2)) ' out of ' num2str(nfeats) ' features (diversity).']);
else
    out.pctage = ones(1,nfeats);
    out.selvars = true(1,nfeats);
end

end