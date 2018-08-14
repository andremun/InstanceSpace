function [X,Y,out] = autoNormalize(X,Y)

X = bsxfun(@minus,X,min(X,[],1))+1;
out.lambdaX = zeros(1,size(X,2));
for i=1:nfeats
    [X(:,i), out.lambdaX(i)] = boxcox(X(:,i));
end
[X, out.muX, out.sigmaX] = zscore(X);

out.lambdaY = zeros(1,size(Y,2));
for i=1:nalgos
    [Y(:,i), out.lambdaY(i)] = boxcox(Y(:,i));
end
[Y, out.muY, out.sigmaY] = zscore(Y);
end