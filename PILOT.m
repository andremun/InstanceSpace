function out = PILOT(X, Y, featlabels, opts)
% -------------------------------------------------------------------------
% PILOT.m
% -------------------------------------------------------------------------
%
% By: Mario Andres Munoz Acosta
%     School of Mathematics and Statistics
%     The University of Melbourne
%     Australia
%     2020
%
% -------------------------------------------------------------------------

errorfcn = @(alpha,Xbar,n,m) nanmean(nanmean((Xbar-(reshape(alpha((2*n)+1:end),m,2)*... % B,C
                                                    reshape(alpha(1:2*n),2,n)...        % A
                                                   *Xbar(:,1:n)')').^2,1),2);
n = size(X, 2); % Number of features
Xbar = [X Y];
m = size(Xbar, 2);
Hd = pdist(X)';
if exist('gcp','file')==2
    mypool = gcp('nocreate');
    if ~isempty(mypool)
        nworkers = mypool.NumWorkers;
    else
        nworkers = 0;
    end
else
    nworkers = 0;
end

if opts.analytic
    disp('  -> PILOT is solving analyticaly the projection problem.');
    disp('  -> This won''t take long.');
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
    if isfield(opts,'alpha') && isnumeric(opts.alpha) && ...
                size(opts.alpha,1)==2*m+2*n && size(opts.alpha,2)==1
        disp('  -> PILOT is using a pre-calculated solution.');
        idx = 1;
        out.alpha = opts.alpha;
    else
        if isfield(opts,'X0') && isnumeric(opts.X0) && ...
                size(opts.X0,1)==2*m+2*n && size(opts.X0,2)>=1
            disp('  -> PILOT is using a user defined starting points for BFGS.');
            X0 = opts.X0;
            opts.ntries = size(opts.X0,2);
        else
            disp('  -> PILOT is using a random starting points for BFGS.');
            state = rng;
            rng('default');
            X0 = 2*rand(2*m+2*n, opts.ntries)-1;
            rng(state);
        end
        alpha = zeros(2*m+2*n, opts.ntries);
        eoptim = zeros(1, opts.ntries);
        perf = zeros(1, opts.ntries);
        disp('-------------------------------------------------------------------------');
        disp('  -> PILOT is solving numerically the projection problem.');
        disp('  -> This may take a while. Trials will not be run sequentially.');
        disp('-------------------------------------------------------------------------');
        parfor (i=1:opts.ntries,nworkers)
            [alpha(:,i),eoptim(i)] = fminunc(errorfcn, X0(:,i), ...
                                             optimoptions('fminunc','Algorithm','quasi-newton',...
                                                                    'Display','off',...
                                                                    'UseParallel',false),...
                                             Xbar, n, m);
            aux = alpha(:,i);
            A = reshape(aux(1:2*n),2,n);
            Z = X*A';
            perf(i) = corr(Hd,pdist(Z)');
            disp(['    -> PILOT has completed trial ' num2str(i)]);
        end
        out.X0 = X0;
        out.alpha = alpha;
        out.eoptim = eoptim;
        out.perf = perf;
        [~,idx] = max(out.perf);
    end
    out.A = reshape(out.alpha(1:2*n,idx),2,n);
    out.Z = X*out.A';
    B = reshape(out.alpha((2*n)+1:end,idx),m,2);
    Xhat = out.Z*B';
    out.C = B(n+1:m,:)';
    out.B = B(1:n,:);
    out.error = sum(sum((Xbar-Xhat).^2,2));
    out.R2 = diag(corr(Xbar,Xhat)).^2;
end

disp('-------------------------------------------------------------------------');
disp('  -> PILOT has completed. The projection matrix A is:');
out.summary = cell(3, n+1);
out.summary(1,2:end) = featlabels;
out.summary(2:end,1) = {'Z_{1}','Z_{2}'};
out.summary(2:end,2:end) = num2cell(round(out.A,4));
disp(' ');
disp(out.summary);

end