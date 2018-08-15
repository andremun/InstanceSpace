function [X,Y,out] = autoNormalize(X,Y)

nfeats = size(X,2);
X = bsxfun(@minus,X,min(X,[],1))+1;
out.lambdaX = zeros(1,nfeats);
for i=1:nfeats
    [X(:,i), out.lambdaX(i)] = boxcox(X(:,i));
end
[X, out.muX, out.sigmaX] = zscore(X);

nalgos = size(Y,2);
out.lambdaY = zeros(1,nalgos);
for i=1:nalgos
    [Y(:,i), out.lambdaY(i)] = boxcox(Y(:,i));
end
[Y, out.muY, out.sigmaY] = zscore(Y);
end