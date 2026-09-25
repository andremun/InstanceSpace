function out = CLOISTER(X, A, opts)
% CLOISTER  Estimate the empirical boundary of the instance space.
%
%   out = CLOISTER(X, A, opts)
%
%   Enumerates every combination of feature lower/upper bounds (a
%   hypercube's corners), discards combinations that contradict the
%   significant pairwise feature correlations (opts.corrThreshold,
%   opts.pval), and projects the surviving corners through the PILOT
%   projection matrix A to trace a boundary polygon in the instance
%   space. If the feature count exceeds opts.maxFeatures, the corner
%   enumeration is skipped and a plain convex hull of the projected
%   instances is used instead.
%
%   Inputs
%     X    - (ninst x nfeats) feature matrix; may contain sparse NaNs
%     A    - (ndim x nfeats) PILOT projection matrix (model.pilot.A)
%     opts - struct with fields (see ISAdefaults for defaults):
%              pval          double  significance level for feature correlations (0.05)
%              corrThreshold double  |correlation| above which a corner combination
%                                    contradicting the trend is discarded (0.70)
%              maxFeatures   int     feature-count guard before falling back
%                                    to a plain convex hull (20)
%
%   Outputs
%     out  - struct with fields:
%              Zedge   boundary vertices using every corner. 2D: a closed
%                      polygon (first vertex repeated last). 3D: the
%                      convex hull's vertices
%              Zecorr  boundary vertices using only correlation-consistent
%                      corners (same as Zedge when the correlation
%                      threshold rejects nothing, or when maxFeatures
%                      triggers the convex-hull fallback)
%              ZedgeFaces, ZecorrFaces
%                      3D only: (nfaces x 3) hull triangulation, indices
%                      into the rows of Zedge/Zecorr. Empty in 2D

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
% -------------------------------------------------------------------------

fprintf('[CLOISTER] CLOISTER is using correlation to estimate a boundary for the space.\n');

nfeats = size(X,2);
[rho,pval] = corr(X);
rho = rho.*(pval<opts.pval);

% omitnan: a feature column can still carry sparse NaNs here (buildIS only
% drops instances where every feature is NaN, and features whose NaN
% fraction exceeds opts.prelim.nanThreshold -- a column under that
% threshold keeps its remaining NaN entries through PRELIM/SIFTED/PILOT).
% Without omitnan, min/max would return NaN for that column and propagate
% through Xedge/Zedge into convhull, which errors on NaN input.
Xbnds = [min(X,[],1,'omitnan'); max(X,[],1,'omitnan')];
% Guard: if too many features, the bit-matrix enumeration below would
% produce an intractable matrix. Use convex hull of Z as a safe fallback.
if ~isfield(opts, 'maxFeatures'), opts.maxFeatures = 20; end
MAX_FEATS = opts.maxFeatures;
if nfeats > MAX_FEATS
    warning('ISA:CLOISTER:tooManyFeatures', ...
        'CLOISTER skipped: %d features exceeds limit of %d. Using convex hull as boundary.', ...
        nfeats, MAX_FEATS);
    [out.Zedge, out.ZedgeFaces] = boundaryHull(X*A');
    out.Zecorr      = out.Zedge;
    out.ZecorrFaces = out.ZedgeFaces;
    fprintf('[CLOISTER] CLOISTER has completed.\n');
    return;
end
% Pure-MATLAB replacement for de2bi (no Communications Toolbox required)
idx = rem(floor((0:2^nfeats-1)' .* 2.^(-(nfeats-1:-1:0))), 2) + 1;
ncomb = size(idx,1);
Xedge = zeros(ncomb,nfeats);
remove = false(ncomb,1);
for i=1:ncomb
   ind = sub2ind([2 nfeats],idx(i,:),1:nfeats);
   Xedge(i,:) = Xbnds(ind)';
   for j=1:nfeats
       for k=j+1:nfeats
           % Check for valid points give the correlation trend
           if rho(j,k)>opts.corrThreshold && sign(Xedge(i,j))~=sign(Xedge(i,k))
               remove(i) = true;
           elseif rho(j,k)<-opts.corrThreshold && sign(Xedge(i,j))==sign(Xedge(i,k))
               remove(i) = true;
           end
           if remove(i)
               break;
           end
       end
       if remove(i)
           break;
       end
   end
end
[out.Zedge, out.ZedgeFaces] = boundaryHull(Xedge*A');

try
    [out.Zecorr, out.ZecorrFaces] = boundaryHull(Xedge(~remove,:)*A');
catch
    fprintf('[CLOISTER] The acceptable correlation threshold was too strict.\n');
    fprintf('[CLOISTER] The features are weakly correlated.\n');
    fprintf('[CLOISTER] Please consider increasing it.\n');
    out.Zecorr      = out.Zedge;
    out.ZecorrFaces = out.ZedgeFaces;
end
fprintf('[CLOISTER] CLOISTER has completed.\n');
end

% =========================================================================
function [V, F] = boundaryHull(Z)
% Convex hull of the projected points in the projection's own
% dimensionality (#50). Before, the hull always used Z(:,1:2), so a 3D
% projection got a boundary that ignored its third coordinate.
%   2D: V is the closed polygon (first vertex repeated last, as convhull
%       returns it) and F is empty -- the same Zedge as before.
%   3D: V holds the hull's unique vertices and F is the (nfaces x 3)
%       triangulation, with indices into the rows of V. Coplanar points
%       give a flat, fan-triangulated polygon.
if size(Z, 2) == 3 && rank(Z - mean(Z, 1)) < 3
    % Coplanar points (e.g. a 3D projection of only two features) have no
    % volumetric hull, and convhull would error. Their hull is a flat
    % polygon: find it in the plane's own 2D coordinates and triangulate it
    % as a fan, so the output keeps the 3D vertices-plus-faces format.
    Zc = Z - mean(Z, 1);
    [~, ~, basis] = svd(Zc, 'econ');
    P = Zc*basis(:, 1:2);
    K = convhull(P(:,1), P(:,2));
    K = K(1:end-1);                 % drop the repeated closing vertex
    V = Z(K, :);
    n = numel(K);
    F = [ones(n-2, 1), (2:n-1)', (3:n)'];
elseif size(Z, 2) == 3
    K = convhull(Z(:,1), Z(:,2), Z(:,3));
    [vidx, ~, F] = unique(K(:));
    V = Z(vidx, :);
    F = reshape(F, size(K));
else
    K = convhull(Z(:,1), Z(:,2));
    V = Z(K, :);
    F = [];
end
end
