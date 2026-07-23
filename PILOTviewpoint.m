function out = PILOTviewpoint(Z, Y, opts)
% PILOTviewpoint  Find the optimal 2D camera viewpoint(s) of a 3D PILOT
% projection (spec §5.2; Equation 2 of Simpson et al., 2025).
%
%   out = PILOTviewpoint(Z, Y, opts)
%
%   Inputs
%     Z    - (ninst x 3) 3D PILOT projection (opts.pilot.dims == 3)
%     Y    - (ninst x nalgos) performance matrix
%     opts - struct with fields:
%              viewGroups  cell array of algorithm index vectors, one
%                          viewpoint computed per group (default {} ->
%                          one global viewpoint over all algorithms)
%
%   Outputs
%     out  - struct with fields:
%              groups     the resolved cell array of algorithm groups
%              V          (ngroups x 1) cell array of fitted 3x2 view
%                         matrices [v1 v2]
%              azimuth    (ngroups x 1) azimuth angle (radians) of the
%                         viewing direction cross(v1,v2), for view(az,el)
%              elevation  (ngroups x 1) elevation angle (radians)
%
%   For each group, finds V = [v1 v2] (3x2) minimising
%       ||Y_group - (Z*V)*M||_F^2 + LAMBDA*dot(v1,v2)
%   where M is the OLS regression of Y_group onto the flattened 2D view
%   Z*V (concentrated out analytically for each candidate V, since
%   Equation 2 only minimises over V). The dot-product penalty encourages
%   v1 and v2 towards orthogonality; LAMBDA=0.2 is the paper-calibrated
%   weight (not user-exposed, per spec §5.2). After optimisation, v1 and
%   v2 are each rescaled to unit magnitude once, and the viewing
%   direction is their cross product, converted to azimuth/elevation.
%
%   By: Mario Andres Munoz Acosta
%       School of Mathematics and Statistics
%       The University of Melbourne
%       Australia
%       2026

narginchk(3, 3);
if size(Z,2) ~= 3
    error('ISA:PILOTviewpoint:not3D', ...
        'Z must be a (ninst x 3) 3D PILOT projection.');
end
if size(Z,1) ~= size(Y,1)
    error('ISA:PILOTviewpoint:sizeMismatch', ...
        'Z and Y must have the same number of rows (got %d vs %d).', ...
        size(Z,1), size(Y,1));
end
if ~isfield(opts, 'viewGroups') || isempty(opts.viewGroups)
    groups = {1:size(Y,2)};
else
    groups = opts.viewGroups;
end

LAMBDA = 0.2; % paper-calibrated orthogonality penalty weight (spec 5.2)

ngroups = numel(groups);
out.groups    = groups;
out.V         = cell(ngroups, 1);
out.azimuth   = zeros(ngroups, 1);
out.elevation = zeros(ngroups, 1);

for g = 1:ngroups
    Yg = Y(:, groups{g});
    V0 = [1 0; 0 1; 0 0]; % start from the Z1-Z2 plane
    theta = fminunc(@(theta) pilotViewErrorFcn(theta, Z, Yg, LAMBDA), V0(:), ...
                     optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'off'));
    V = reshape(theta, 3, 2);
    v1 = V(:,1) ./ norm(V(:,1));
    v2 = V(:,2) ./ norm(V(:,2));
    out.V{g} = [v1 v2];
    viewdir = cross(v1, v2);
    [out.azimuth(g), out.elevation(g)] = cart2sph(viewdir(1), viewdir(2), viewdir(3));
end

end
% =========================================================================
function err = pilotViewErrorFcn(theta, Z, Y, lambda)
% Cost for a candidate view matrix V = reshape(theta,3,2): the flattened
% 2D view Z*V is regressed onto Y via OLS (concentrated out, since
% Equation 2 minimises over V only), and the orthogonality penalty is
% applied to the raw (unnormalised) columns of V.
V = reshape(theta, 3, 2);
Zview = Z * V;         % (ninst x 2)
M = Zview \ Y;          % (2 x nalgos) OLS coefficients
Yhat = Zview * M;        % (ninst x nalgos)
err = nanmean(nanmean((Y - Yhat).^2, 1), 2) + lambda * dot(V(:,1), V(:,2));
end
