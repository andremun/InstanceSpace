function out = PBLDR(X, Y, opts)

errorfcn = @(alpha,X,n,m) mean(mean((X-(reshape(alpha((2*n)+1:end),m,2)*... % B,C
                                        reshape(alpha(1:2*n),2,n)...        % A
                                        *X(:,1:n)')').^2,1),2);

n = size(X, 2); % Number of features
Xbar = [X Y];
m = size(Xbar, 2);
Hd = pdist(X)';

if opts.analytic
    Xbar = Xbar';
    X = X';
    [V,D] = eig(Xbar*Xbar');
    [~,idx] = sort(abs(diag(D)),'descend');
    V = V(:,idx(1:2));
    out.B = V(1:n,:);
    out.C = V(n+1:m,:)';
    Xr = X'/(X*X');
    out.A = V'*Xbar*Xr;
    out.Z = out.A*X;
    Xhat = [out.B*out.Z; out.C'*out.Z];
    out.error = sum(sum((Xbar-Xhat).^2,2));
    out.R2 = diag(corr(Xbar',Xhat')).^2;
else
    alpha = zeros(2*m+2*n, opts.ntries);
    out.eoptim = zeros(1, opts.ntries);
    out.perf = zeros(1, opts.ntries);
    initstr = ['2*rand(' num2str(2*m+2*n) ',1)-1'];

    for i=1:opts.ntries
        [alpha(:,i),out.eoptim(i)] = bipopcmaes(errorfcn, ...
                                                initstr, ...
                                                1, ...
                                                opts.cmaopts, ...
                                                Xbar, ...
                                                n, ...
                                                m);
        A = reshape(alpha(1:2*n,i),2,n);
        Z = X*A';
        out.perf(i) = corr(Hd,pdist(Z)');
    end

    [~,idx] = max(out.perf);
    out.A = reshape(alpha(1:2*n,idx),2,n);
    out.Z = X*out.A';
    B = reshape(alpha((2*n)+1:end,idx),m,2);
    Xhat = out.Z*B';
    out.C = B(n+1:m,:)';
    out.B = B(1:n,:);
    out.error = sum(sum((Xbar-Xhat).^2,2));
    out.R2 = diag(corr(Xbar,Xhat)).^2;
end
end