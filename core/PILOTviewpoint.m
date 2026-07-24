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
function out = PILOTviewpoint(Z, Y, opts)
% PILOTviewpoint  Find the optimal 2D camera viewpoint(s) of a 3D PILOT
% projection (spec §5.2; Equation 2 of Simpson et al., 2025).
%
%   out = PILOTviewpoint(Z, Y, opts)
%
%   Inputs
%     Z    - (ninst x n) 3D PILOT projection, n==3 (opts.pilot.dims == 3)
%     Y    - (ninst x nalgos) performance matrix
%     opts - struct with fields:
%              viewGroups  cell array of algorithm index vectors, one
%                          viewpoint computed per group (default {} ->
%                          one global viewpoint over all algorithms)
%              ntries      BFGS multi-start restarts (default 10, same
%                          convention as opts.pilot.ntries)
%              X0          optional (2n+2*numel(group)) x ntries
%                          user-supplied starting points, reused for
%                          every group whose size matches
%
%   Outputs
%     out  - struct with fields:
%              groups     the resolved cell array of algorithm groups
%              A          (ngroups x 1) cell array of fitted 2xn view
%                         matrices [row1=v1; row2=v2] flattening Z onto
%                         the viewing plane
%              azimuth    (ngroups x 1) azimuth angle (radians) of the
%                         viewing direction cross(v1,v2), for view(az,el)
%              elevation  (ngroups x 1) elevation angle (radians)
%
%   For each group, jointly fits A (2xn, the view) and C (nalgos x 2, the
%   performance reconstruction) via BFGS, minimising
%       ||Y_group - (C*A*Z')'||_F^2 + LAMBDA*|dot(v1,v2)|
%   where v1=A(1,:), v2=A(2,:), following the same multi-start /
%   topological-preservation trial-selection scheme as PILOT.m's
%   numerical branch (Hd = pdist(Z), best trial = highest corr(Hd,
%   pdist(Z*A'))). LAMBDA=0.2 is the paper-calibrated orthogonality
%   penalty weight (not user-exposed, per spec §5.2). v1 and v2 are
%   rescaled to unit magnitude once per trial (both for the
%   topological-preservation scoring and the stored solution), replacing
%   the original PILOTANGLE.m reference's redundant and axis-wrong
%   sum(A.^2') row-normalisation applied twice.

narginchk(3, 3);
if size(Z,1) ~= size(Y,1)
    error('ISA:PILOTviewpoint:sizeMismatch', ...
        'Z and Y must have the same number of rows (got %d vs %d).', ...
        size(Z,1), size(Y,1));
end
if size(Z,2) ~= 3
    % cross(v1,v2)/cart2sph below are only meaningful for a 3D viewing
    % direction; fail fast here instead of an index-out-of-bounds on
    % viewdir(3) once a trial has already been fitted.
    error('ISA:PILOTviewpoint:not3D', ...
        'Z must be a (ninst x 3) 3D PILOT projection (got %d columns).', size(Z,2));
end
if ~isfield(opts, 'ntries'), opts.ntries = 10; end
if ~isfield(opts, 'viewGroups') || isempty(opts.viewGroups)
    groups = {1:size(Y,2)};
elseif iscell(opts.viewGroups)
    groups = opts.viewGroups;
else
    % jsondecode collapses a JSON array of equal-length arrays into a
    % plain numeric matrix instead of a cell array, so opts.viewGroups
    % may arrive this way after an options.json round-trip. Treat each
    % row as one group.
    groups = num2cell(opts.viewGroups, 2);
end

LAMBDA = 0.2; % paper-calibrated orthogonality penalty weight (spec 5.2)

n = size(Z, 2); % projection dimensionality being viewed (3 for 3D PILOT)
Hd = pdist(Z)';
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

ngroups = numel(groups);
out.groups    = groups;
out.A         = cell(ngroups, 1);
out.azimuth   = zeros(ngroups, 1);
out.elevation = zeros(ngroups, 1);

for g = 1:ngroups
    Yg = Y(:, groups{g});
    n2 = size(Yg, 2);
    errorfcn = @(theta) pilotViewErrorFcn(theta, Z, Yg, n, n2, LAMBDA);

    if isfield(opts,'X0') && isnumeric(opts.X0) && ...
            size(opts.X0,1)==2*n+2*n2 && size(opts.X0,2)>=1
        X0 = opts.X0;
        ntries = size(opts.X0,2);
    else
        state = rng;
        rng('default');
        X0 = 2*rand(2*n+2*n2, opts.ntries)-1;
        rng(state);
        ntries = opts.ntries;
    end

    theta = zeros(2*n+2*n2, ntries);
    perf  = zeros(1, ntries);
    parfor (i=1:ntries, nworkers)
        thetai = fminunc(errorfcn, X0(:,i), ...
                          optimoptions('fminunc','Algorithm','quasi-newton',...
                                                 'Display','off',...
                                                 'UseParallel',false,...
                                                 'MaxIterations',30000,...
                                                 'FunctionTolerance',1e-20));
        A = reshape(thetai(1:2*n), 2, n);
        A(1,:) = A(1,:) / max(norm(A(1,:)), eps);
        A(2,:) = A(2,:) / max(norm(A(2,:)), eps);
        thetai(1:2*n) = A(:);
        theta(:,i) = thetai;
        perf(i) = corr(Hd, pdist(Z*A')');
    end
    [~,idx] = max(perf);
    A = reshape(theta(1:2*n,idx), 2, n);
    out.A{g} = A;
    viewdir = cross(A(1,:), A(2,:));
    [out.azimuth(g), out.elevation(g)] = cart2sph(viewdir(1), viewdir(2), viewdir(3));
end

end
% =========================================================================
function err = pilotViewErrorFcn(theta, Z, Y, n, n2, lambda)
% Joint view (A, 2xn) + performance-reconstruction (C, n2x2) cost:
%   ||Y - (C*A*Z')'||_F^2 + lambda*|dot(v1,v2)|
% where v1,v2 are A's rows RESCALED TO UNIT NORM. The reconstruction term
% C*A*Z' is invariant to rescaling A's rows while inversely rescaling C's
% columns, so penalising the raw (unnormalised) dot product -- as
% PILOTANGLE.m's reference does -- lets the optimiser shrink the penalty
% toward 0 just by shrinking A, without the two directions actually
% becoming orthogonal; normalising only inside the penalty (not the
% reconstruction term) keeps it a true orthogonality measure regardless
% of how BFGS happens to split scale between A and C.
A = reshape(theta(1:2*n), 2, n);
C = reshape(theta((2*n)+1:end), n2, 2);
Yhat = (C*A*Z')';
v1 = A(1,:) / max(norm(A(1,:)), eps);
v2 = A(2,:) / max(norm(A(2,:)), eps);
err = nanmean(nanmean((Y-Yhat).^2,1),2) + lambda*abs(dot(v1,v2));
end
