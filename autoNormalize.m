function [X,Y,out] = autoNormalize(X, Y)

nfeats = size(X,2);
nalgos = size(Y,2);
disp('-> Auto-normalizing the data using Box-Cox and Z transformations.');
out.minX = min(X,[],1);
X = bsxfun(@minus,X,out.minX)+1;
out.lambdaX = zeros(1,nfeats);
out.muX = zeros(1,nfeats);
out.sigmaX = zeros(1,nfeats);
for i=1:nfeats
    aux = X(:,i);
    idx = isnan(aux);
    [aux, out.lambdaX(i)] = boxcox(aux(~idx));
    [aux, out.muX(i), out.sigmaX(i)] = zscore(aux);
    X(~idx,i) = aux;
end

out.minY = min(Y(:));
Y = (Y-out.minY)+eps;
out.lambdaY = zeros(1,nalgos);
out.muY = zeros(1,nalgos);
out.sigmaY = zeros(1,nalgos);
for i=1:nalgos
    aux = Y(:,i);
    idx = isnan(aux);
    [aux, out.lambdaY(i)] = boxcox(aux(~idx));
    [aux, out.muY(i), out.sigmaY(i)] = zscore(aux);
    Y(~idx,i) = aux;
end

end