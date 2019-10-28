function out = PBLDR(X, Y, opts)

errorfcn = @(alpha,X,n,m) mean(mean((X-(reshape(alpha((2*n)+1:end),m,2)*... % B,C
                                        reshape(alpha(1:2*n),2,n)...        % A
                                        *X(:,1:n)')').^2,1),2);

% Here there is a problem with Inf values. It won't be able to calculate
% the projection. You need to do something different.
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
    out.alpha = zeros(2*m+2*n, opts.ntries);
    out.eoptim = zeros(1, opts.ntries);
    out.perf = zeros(1, opts.ntries);
    out.X0 = out.alpha;
    for i=1:opts.ntries
        out.X0(:,i) = 2*rand(2*m+2*n,1)-1;
        [out.alpha(:,i),out.eoptim(i)] = fminunc(errorfcn, out.X0(:,i), ...
                                                 optimoptions('fminunc','Algorithm','quasi-newton',...
                                                                        'Display','off',...
                                                                        'UseParallel',false),...
                                                 Xbar, n, m);
        A = reshape(out.alpha(1:2*n,i),2,n);
        Z = X*A';
        out.perf(i) = corr(Hd,pdist(Z)');
        disp(['    -> Trial number: ' num2str(i) ' completed']);
    end

    [~,idx] = max(out.perf);
    out.A = reshape(out.alpha(1:2*n,idx),2,n);
    out.Z = X*out.A';
    B = reshape(out.alpha((2*n)+1:end,idx),m,2);
    Xhat = out.Z*B';
    out.C = B(n+1:m,:)';
    out.B = B(1:n,:);
    out.error = sum(sum((Xbar-Xhat).^2,2));
    out.R2 = diag(corr(Xbar,Xhat)).^2;
end
end