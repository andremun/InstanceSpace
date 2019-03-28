function [X,Y,out] = autoNormalize(X, Y, opts)

nfeats = size(X,2);
nalgos = size(Y,2);
if opts.flag
    out.minX = min(X,[],1);
    X = bsxfun(@minus,X,out.minX)+1;
    out.lambdaX = zeros(1,nfeats);
    for i=1:nfeats
        [X(:,i), out.lambdaX(i)] = boxcox(X(:,i));
    end
    [X, out.muX, out.sigmaX] = zscore(X);
    
    Y(Y==0) = eps;
    out.lambdaY = zeros(1,nalgos);
    for i=1:nalgos
        [Y(:,i), out.lambdaY(i)] = boxcox(Y(:,i));
    end
    [Y, out.muY, out.sigmaY] = zscore(Y);
else
    out.lambdaX = ones(1,nfeats);
    out.muX = zeros(1,nfeats);
    out.sigmaX = ones(1,nfeats);
    out.lambdaY = ones(1,nalgos);
    out.muY = zeros(1,nalgos);
    out.sigmaY = ones(1,nalgos);
end

end