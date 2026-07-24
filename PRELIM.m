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
%   Simpson, C., Munoz, M.A., Kandanaarachchi, S. & Campello, R.J.G.B.
%   (2025). ISA3: A 3-dimensional expansion of Instance Space Analysis.
%   Machine Learning, 114, 240. https://doi.org/10.1007/s10994-025-06871-5
%
%   Smith-Miles, K. & Munoz, M.A. (2023). Instance Space Analysis for
%   Algorithm Testing. ACM Computing Surveys, 55(12), Article 255.
%   https://doi.org/10.1145/3572895
% -------------------------------------------------------------------------
function [X, Y, out] = PRELIM(X, Y, opts)
% PRELIM  Pre-process instance-space data before projection.
%
%   [X, Y, out] = PRELIM(X, Y, opts)
%
%   Inputs
%     X     - (ninst x nfeats) feature matrix; may contain NaN
%     Y     - (ninst x nalgos) performance matrix; may contain NaN
%     opts  - struct with fields:
%               MaxPerf       logical  true = maximise performance (default false)
%               AbsPerf       logical  true = absolute threshold (default false)
%               epsilon       double   good-performance threshold (default 0.05)
%               betaThreshold double   easy-instance fraction (default 0.55)
%               auto          logical  run auto pre-processing (default true)
%               bound         logical  bound outliers (default true)
%               norm          logical  Box-Cox + Z normalisation (default true)
%               iqrMultiplier double   outlier bound = median +/- N*IQR (default 5)
%
%   Outputs
%     X     - pre-processed feature matrix
%     Y     - pre-processed performance matrix
%     out   - struct with fields:
%               Ybest, Ybin, P, numGoodAlgos, beta  (performance summary)
%               medval, iqrange, hibound, lobound    (outlier bounds)
%               minX, lambdaX, muX, sigmaX           (feature normalisation)
%               minY, lambdaY, muY, sigmaY           (performance normalisation)

narginchk(3, 3);
if size(X, 1) ~= size(Y, 1)
    error('ISA:PRELIM:sizeMismatch', ...
        'X and Y must have the same number of rows (got %d vs %d).', ...
        size(X,1), size(Y,1));
end
if ~isstruct(opts)
    error('ISA:PRELIM:badOpts', 'opts must be a struct.');
end
if ~isfield(opts, 'iqrMultiplier'), opts.iqrMultiplier = 5; end

Yraw = Y;
nalgos = size(Y, 2);
% -------------------------------------------------------------------------
% Determine whether the performance of an algorithm is a cost measure to
% be minimized or a profit measure to be maximized. Moreover, determine
% whether we are using an absolute threshold as good performance (the
% algorithm has a performance better than the threshold) or a relative
% performance (the algorithm has a performance that is similar to the
% best algorithm minus a percentage).
fprintf('[PRELIM] Calculating the binary measure of performance\n');
msg = 'An algorithm is good if its performance is ';
if opts.MaxPerf
    Yaux = Y;
    Yaux(isnan(Yaux)) = -Inf;
    [out.Ybest, out.P] = max(Yaux, [], 2);
    YbestTie = out.Ybest;  % pre-eps-substitution snapshot, used for tie detection below
    if opts.AbsPerf
        out.Ybin = Yaux >= opts.epsilon;
        msg = [msg 'higher than ' num2str(opts.epsilon)];
    else
        if mean(out.Ybest==0) > 0.05
            warning('ISA:PRELIM:manyZeroBest', ...
                ['More than 5%% of instances have a best-algorithm performance of ' ...
                 'exactly zero; the relative-performance matrix will be close to 1 ' ...
                 'everywhere for these instances.']);
        end
        out.Ybest(out.Ybest==0) = eps;
        Y(Y==0) = eps;
        Y = 1 - bsxfun(@rdivide, Y, out.Ybest);
        out.Ybin = (1 - bsxfun(@rdivide, Yaux, out.Ybest)) <= opts.epsilon;
        msg = [msg 'within ' num2str(round(100.*opts.epsilon)) '% of the best.'];
    end
else
    Yaux = Y;
    Yaux(isnan(Yaux)) = Inf;
    [out.Ybest, out.P] = min(Yaux, [], 2);
    YbestTie = out.Ybest;  % pre-eps-substitution snapshot, used for tie detection below
    if opts.AbsPerf
        out.Ybin = Yaux <= opts.epsilon;
        msg = [msg 'less than ' num2str(opts.epsilon)];
    else
        if mean(out.Ybest==0) > 0.05
            warning('ISA:PRELIM:manyZeroBest', ...
                ['More than 5%% of instances have a best-algorithm performance of ' ...
                 'exactly zero; the relative-performance matrix will be close to 1 ' ...
                 'everywhere for these instances.']);
        end
        out.Ybest(out.Ybest==0) = eps;
        Y(Y==0) = eps;
        Y = bsxfun(@rdivide, Y, out.Ybest) - 1;
        out.Ybin = (bsxfun(@rdivide, Yaux, out.Ybest) - 1) <= opts.epsilon;
        msg = [msg 'within ' num2str(round(100.*opts.epsilon)) '% of the best.'];
    end
end
fprintf('[PRELIM] %s\n', msg);
% -------------------------------------------------------------------------
% Testing for ties. If there is a tie in performance, we pick an algorithm
% at random. Compared against YbestTie (captured before the eps
% substitution above) so exact-zero best scores are still matched correctly.
bestAlgos = bsxfun(@eq, Yraw, YbestTie);
multipleBestAlgos = sum(bestAlgos, 2) > 1;
aidx = 1:nalgos;
tieRows = find(multipleBestAlgos)';  % typically very few rows
for i = tieRows
    aux = aidx(bestAlgos(i,:));
    out.P(i) = aux(randi(numel(aux)));
end
fprintf('[PRELIM] For %s%% of the instances there is more than one best algorithm. Random selection is used to break ties.\n', ...
    num2str(round(100.*mean(multipleBestAlgos))));
out.numGoodAlgos = sum(out.Ybin, 2);
out.beta = out.numGoodAlgos > (opts.betaThreshold * nalgos);
% -------------------------------------------------------------------------
if opts.auto
    fprintf('[PRELIM] Auto-pre-processing.\n');
end
out.medval  = nanmedian(X, 1);
out.iqrange = iqr(X, 1);
out.hibound = out.medval + opts.iqrMultiplier.*out.iqrange;
out.lobound = out.medval - opts.iqrMultiplier.*out.iqrange;
if opts.auto && opts.bound
    fprintf('[PRELIM] Removing extreme outliers from the feature values.\n');
    himask = bsxfun(@gt, X, out.hibound);
    lomask = bsxfun(@lt, X, out.lobound);
    X = X.*~(himask | lomask) + bsxfun(@times, himask, out.hibound) + ...
                                bsxfun(@times, lomask, out.lobound);
end

nfeats = size(X, 2);
nalgos = size(Y, 2);
out.minX    = min(X, [], 1, 'omitnan');
out.lambdaX = zeros(1, nfeats);
out.muX     = zeros(1, nfeats);
out.sigmaX  = zeros(1, nfeats);
out.minY    = nanmin(Y(:));
out.lambdaY = zeros(1, nalgos);
out.muY     = zeros(1, nalgos);
out.sigmaY  = zeros(1, nalgos);
if opts.auto && opts.norm
    fprintf('[PRELIM] Auto-normalizing the data using Box-Cox and Z transformations.\n');
    X = bsxfun(@minus, X, out.minX) + 1;
    for i = 1:nfeats
        aux = X(:, i);
        idx = isnan(aux);
        [aux, out.lambdaX(i)] = boxcox(aux(~idx));
        [aux, out.muX(i), out.sigmaX(i)] = zscore(aux);
        X(~idx, i) = aux;
    end

    Y = (Y - out.minY) + eps;
    for i = 1:nalgos
        aux = Y(:, i);
        idx = isnan(aux);
        [aux, out.lambdaY(i)] = boxcox(aux(~idx));
        [aux, out.muY(i), out.sigmaY(i)] = zscore(aux);
        Y(~idx, i) = aux;
    end
end

end
