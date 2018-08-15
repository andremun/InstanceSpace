function [X, out] = checkDiversity(X, DIVTHRESHOLD)

[ninst,nfeats] = size(X);
out.pctage = zeros(1,nfeats);
for i=1:nfeats
    out.pctage(i) = length(unique(X(:,i)))./ninst;
end
out.selvars = out.pctage>DIVTHRESHOLD;
X = X(:,out.selvars);
end