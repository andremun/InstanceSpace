function [X, out] = checkDiversity(X, DIVTHRESHOLD)

nfeats = size(X,2);
out.pctage = zeros(1,nfeats);
for i=1:nfeats
    out.pctage(i) = length(unique(X(:,i)))./ninst;
end
out.selvars = out.pctage>DIVTHRESHOLD;
X = X(:,selvars);
end