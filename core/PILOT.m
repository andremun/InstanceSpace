function out = PILOT(X, Y, featlabels, opts)
% PILOT  Project features onto a 2D or 3D instance space (Munoz et al.,
% Mach Learn 2018), finding A/B/C such that Z=X*A' and [X Y] is
% reconstructed from Z as closely as possible.
%
%   out = PILOT(X, Y, featlabels, opts)
%
%   Inputs
%     X          - (ninst x nfeats) feature matrix
%     Y          - (ninst x nalgos) performance matrix
%     featlabels - (1 x nfeats) cell array of feature name strings
%     opts       - struct with fields (see ISAdefaults for defaults):
%                    dims      int    projection dimensionality, 2 or 3
%                    method    char   'standard' (BFGS/analytic, below) or
%                                     'pls' (Partial Least Squares via
%                                     plsregress; opts.alpha does not apply)
%                    analytic  logical use the closed-form eigenvector
%                                     solution (method='standard' only);
%                                     falls back to numerical if X is
%                                     rank-deficient
%                    ntries    int    BFGS multi-start restarts (numerical
%                                     branch only)
%                    alpha     double performance-reconstruction cost
%                                     weight (method='standard' only):
%                                     min ||Xtilde-BrZ||^2 +
%                                     alpha*||Y-CrZ||^2
%                    topoWeight double reserved for future use; not wired
%                                     into the cost function
%                    verbose   logical per-trial progress output
%                    precalcAlpha (optional) pre-computed full BFGS
%                                     solution vector, skips optimisation
%                    X0        (optional) user-supplied BFGS starting points
%
%   Outputs
%     out  - struct with fields:
%              A       (dims x nfeats) projection matrix, Z = X*A'
%              B, C    reconstruction matrices for the feature/performance
%                      blocks of [X Y] from Z
%              Z       (ninst x dims) projected instance coordinates
%              error   sum of squared reconstruction error
%              R2      per-column R^2 of the reconstruction
%              summary cell array display of A with feature labels
%              alpha, X0, eoptim, perf  (numerical branch only) raw BFGS
%                      trial results, used to pick the best of opts.ntries

% -------------------------------------------------------------------------
% Instance Space Analysis (ISA) Toolkit
% Copyright (c) 2026 Mario Andres Munoz Acosta and contributors
% School of Computing and Information Systems
% The University of Melbourne, Australia
%
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
% License: https://polyformproject.org/licenses/noncommercial/1.0.0/
%
% You may use, modify, and redistribute this software for non-commercial
% research and educational purposes only. Commercial use requires prior
% written permission. See the LICENSE file for full terms.
%
% Reference:
%   Smith-Miles, K. & Munoz, M.A. (2023). Instance Space Analysis for
%   Algorithm Testing. ACM Computing Surveys, 55(12), Article 255.
%   https://doi.org/10.1145/3572895
%
%   Simpson, C., Munoz, M.A., Kandanaarachchi, S. & Campello, R.J.G.B.
%   (2025). ISA3: A 3-dimensional expansion of Instance Space Analysis.
%   Machine Learning, 114, 240. https://doi.org/10.1007/s10994-025-06871-5
%
%   Munoz, M.A., Villanova, L., Baatar, D. & Smith-Miles, K. (2018).
%   Instance spaces for machine learning classification. Machine
%   Learning, 107(1), 109-147. https://doi.org/10.1007/s10994-017-5629-5
% -------------------------------------------------------------------------

if ~isfield(opts, 'verbose'), opts.verbose = true; end
% Legacy: opts.ISA3D (boolean) -> opts.dims (2|3, spec Appendix A).
if isfield(opts, 'ISA3D') && ~isfield(opts, 'dims')
    opts.dims = 2 + double(logical(opts.ISA3D));
end
if ~isfield(opts, 'dims'),   opts.dims   = 2;          end
if ~isfield(opts, 'alpha'),  opts.alpha  = 1.0;        end % performance-reconstruction cost weight (spec 5.4), 'standard' method only
if ~isfield(opts, 'method'), opts.method = 'standard'; end % 'standard' (BFGS/analytic) or 'pls' (spec 5.3)
if ~any(strcmpi(opts.method, {'standard','pls'}))
    error('ISA:PILOT:invalidMethod', ...
        'opts.method must be ''standard'' or ''pls'' (got ''%s'').', char(string(opts.method)));
end
if ~(isnumeric(opts.dims) && isscalar(opts.dims) && ismember(opts.dims, [2 3]))
    error('ISA:PILOT:invalidDims', ...
        'opts.dims must be 2 or 3 (got %s).', mat2str(opts.dims));
end
% costWeight must stay strictly positive: the analytic branch divides by
% sqrt(costWeight) to unscale C back to true Y units (Inf/NaN at 0,
% complex for negative values), and alpha<=0 has no sensible meaning for
% the numerical branch's weighted loss either.
if ~(isnumeric(opts.alpha) && isscalar(opts.alpha) && isfinite(opts.alpha) && opts.alpha > 0)
    error('ISA:PILOT:invalidAlpha', ...
        'opts.alpha must be a finite, positive scalar (got %s).', mat2str(opts.alpha));
end
d = opts.dims;
costWeight = opts.alpha;

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

if strcmpi(opts.method, 'standard') && opts.analytic && rank(X) < size(X,2)
    warning('ISA:PILOT:rankDeficient', ...
        'Feature matrix rank-deficient; falling back to numerical solution.');
    opts.analytic = false;
end

if strcmpi(opts.method, 'pls')
    % Partial Least Squares alternative (spec §5.3, previously PBLDR).
    % All three projection matrices are natively available from
    % plsregress: the X-weight matrix W gives Ar, the X-loading matrix P
    % gives Br directly, and the Y-loading matrix Q gives Cr. Unlike the
    % standard method, PLS does not require X to be full column rank, and
    % opts.alpha (spec §5.4) does not apply -- plsregress maximises
    % covariance rather than minimising a weighted reconstruction cost.
    %
    % plsregress always mean-centres X and Y internally regardless of
    % whether PRELIM already did so; its scores/loadings are defined in
    % terms of the CENTRED data, so out.error/out.R2 (computed against
    % the true Xbar) must add the means back on reconstruction.
    %
    % out.Z uses plsregress's own XS output directly rather than
    % recomputing (X-Xmean)*stats.W': for the SIMPLS algorithm MATLAB
    % uses, the scores satisfy XS = X0*W directly (SIMPLS deflates the
    % cross-covariance, not X itself, unlike NIPALS -- so no further
    % inv(P'*W) correction applies), but reading XS straight from
    % plsregress guarantees out.Z is consistent with XL/YL by
    % construction regardless of that algorithmic detail.
    if opts.verbose
        fprintf('[PILOT] PILOT is using partial least squares (opts.pilot.method=''pls'').\n');
    end
    Xmean = mean(X, 1);
    Ymean = mean(Y, 1);
    [XL, YL, XS, ~, ~, ~, ~, stats] = plsregress(X, Y, d);
    out.A = stats.W';  % Ar = W' (d x n) -- used by exploreIS to reproject new instances via Z=X*A'
    out.B = XL;          % Br = P (n x d)
    out.C = YL';          % Cr = Q' (d x q)
    out.Z = XS;
    Xhat = out.Z*[out.B; out.C']' + [Xmean Ymean];
    out.error = sum(sum((Xbar-Xhat).^2,2));
    out.R2 = diag(corr(Xbar,Xhat)).^2;
elseif opts.analytic
    % Generalises to 2D or 3D by taking the top d eigenvectors of
    % [Xtilde;Y][Xtilde;Y]' (spec §5.1.1). The performance-reconstruction
    % cost weight (spec §5.4) is folded in via weighted PCA: the Y block is
    % scaled by sqrt(costWeight) before the eigendecomposition, then C is
    % rescaled back so it reconstructs Y in its true (unweighted) units.
    % costWeight=1 reproduces the original unweighted derivation exactly.
    if opts.verbose
        fprintf('[PILOT] PILOT is solving analytically the projection problem.\n');
        fprintf('[PILOT] This won''t take long.\n');
    end
    Xbarw = Xbar;
    Xbarw(:,n+1:m) = sqrt(costWeight) * Xbarw(:,n+1:m);
    XbarwT = Xbarw';          % (m x ninst): weighted features+performance as rows
    XbarT  = Xbar';            % (m x ninst): true (unweighted) features+performance as rows
    Xt     = X';               % (n x ninst): features as rows
    [V,D] = eig(XbarwT*XbarwT');
    [~,idx] = sort(abs(diag(D)),'descend');
    V = V(:,idx(1:d));         % top-d eigenvectors, (m x d)
    out.B = V(1:n,:);           % (n x d)
    out.C = V(n+1:m,:)'./sqrt(costWeight); % (d x q), rescaled back to true Y units
    Xr = Xt'/(Xt*Xt');           % pseudo-inverse of Xt (rank-checked above), (ninst x n)
    out.A = V'*XbarwT*Xr;         % (d x n)
    Zt = out.A*Xt;                 % (d x ninst)
    out.Z = Zt';                   % (ninst x d) -- matches the numerical branch's convention
    Xhat = [out.B*Zt; out.C'*Zt];  % (m x ninst), same orientation as XbarT
    out.error = sum(sum((XbarT-Xhat).^2,2));
    out.R2 = diag(corr(XbarT',Xhat')).^2;
else
    errorfcn = @(theta,Xbar,n,m) pilotErrorFcn(theta,Xbar,n,m,d,costWeight);

    if isfield(opts,'precalcAlpha') && isnumeric(opts.precalcAlpha) && ...
                size(opts.precalcAlpha,1)==d*m+d*n && size(opts.precalcAlpha,2)==1
        if opts.verbose
            fprintf('[PILOT] PILOT is using a pre-calculated solution.\n');
        end
        idx = 1;
        out.alpha = opts.precalcAlpha;
    else
        if isfield(opts,'X0') && isnumeric(opts.X0) && ...
                size(opts.X0,1)==d*m+d*n && size(opts.X0,2)>=1
            if opts.verbose
                fprintf('[PILOT] PILOT is using a user defined starting points for BFGS.\n');
            end
            X0 = opts.X0;
            opts.ntries = size(opts.X0,2);
        else
            if opts.verbose
                fprintf('[PILOT] PILOT is using a random starting points for BFGS.\n');
            end
            state = rng;
            rng('default');
            X0 = 2*rand(d*m+d*n, opts.ntries)-1;
            rng(state);
        end
        alpha = zeros(d*m+d*n, opts.ntries);
        eoptim = zeros(1, opts.ntries);
        perf = zeros(1, opts.ntries);
        if opts.verbose
            fprintf('[PILOT] PILOT is solving numerically the projection problem.\n');
            fprintf('[PILOT] This may take a while. Trials will not be run sequentially.\n');
        end
        parfor (i=1:opts.ntries,nworkers)
            [alphaCol,eoptim(i)] = fminunc(errorfcn, X0(:,i), ...
                                             optimoptions('fminunc','Algorithm','quasi-newton',...
                                                                    'Display','off',...
                                                                    'UseParallel',false),...
                                             Xbar, n, m);
            % alpha(:,i) is a parfor sliced-output variable: it must only
            % ever be written (never read back) with a consistent (:,i)
            % pattern, so A is reshaped from the plain local alphaCol
            % instead of alpha(1:d*n,i) -- that mismatched a :-sliced
            % write against a 1:d*n-ranged read of the same variable,
            % which MATLAB's parfor classifier rejects outright.
            alpha(:,i) = alphaCol;
            A = reshape(alphaCol(1:d*n),d,n);
            Z = X*A';
            perf(i) = corr(Hd,pdist(Z)');
            if opts.verbose
                fprintf('[PILOT] PILOT has completed trial %d\n', i);
            end
        end
        out.X0 = X0;
        out.alpha = alpha;
        out.eoptim = eoptim;
        out.perf = perf;
        [~,idx] = max(out.perf);
    end
    out.A = reshape(out.alpha(1:d*n,idx),d,n);
    out.Z = X*out.A';
    B = reshape(out.alpha((d*n)+1:end,idx),m,d);
    Xhat = out.Z*B';
    out.C = B(n+1:m,:)';
    out.B = B(1:n,:);
    out.error = sum(sum((Xbar-Xhat).^2,2));
    out.R2 = diag(corr(Xbar,Xhat)).^2;
end

out.summary = cell(d+1, n+1);
out.summary(2:end,1) = arrayfun(@(k) sprintf('Z_{%d}',k), 1:d, 'UniformOutput', false);
out.summary(1,2:end) = featlabels;
out.summary(2:end,2:end) = num2cell(round(out.A,4));
if opts.verbose
    fprintf('[PILOT] PILOT has completed. The projection matrix A is:\n\n');
    disp(out.summary);
end

end
% =========================================================================
function err = pilotErrorFcn(theta, Xbar, n, m, d, costWeight)
% Feature+performance reconstruction cost (spec §5.4):
%   ||Ftilde - Br*Z||^2_F + costWeight*||Y - Cr*Z||^2_F
% Xbar columns 1:n are the feature block (Ftilde), n+1:m the performance
% block (Y); costWeight scales only the latter. costWeight=1 reduces to a
% plain nanmean over all m columns, matching the pre-refactor loss exactly.
B = reshape(theta((d*n)+1:end),m,d);
A = reshape(theta(1:d*n),d,n);
Xhat = (B*A*Xbar(:,1:n)')';
sqErr = (Xbar-Xhat).^2;
w = ones(1,m);
w(n+1:m) = costWeight;
err = nanmean(nanmean(bsxfun(@times, sqErr, w), 1), 2);
end
