function [X, out] = checkCorrelation(X,Y,CORTHRESHOLD)

out.rho = corr(X,Y,'rows','pairwise');
[~,row] = sort(abs(out.rho),1,'descend');
out.selvars = unique(row(1:CORTHRESHOLD,:));
X = X(:,out.selvars);
end