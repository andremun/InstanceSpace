function scriptfcn
% scriptfcn  Inject shared CSV-writing/colour-scaling lambdas and plotting
% helpers into the caller's workspace via assignin.
%
%   scriptfcn;
%
%   Called at the top of scriptcsv.m/scriptpng.m/scriptweb.m so they can
%   use writeArray2CSV, writeCell2CSV, makeBndLabels, colorscale,
%   colorscaleg, and the draw*/getPolygonPoints/footprintBoundary/
%   traceAlphaBoundary functions defined below as if they were local
%   functions, without duplicating them in each file. See the individual
%   subfunctions below (axisLimits, applyView, resolveViewAngle,
%   drawSources, drawScatter, drawPortfolioSelections,
%   drawPortfolioFootprint, drawGoodBadFootprint, drawFootprint,
%   drawBinaryPerformance, drawBoundary, getPolygonPoints,
%   footprintBoundary, traceAlphaBoundary) for what each does.

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
% -------------------------------------------------------------------------

writeArray2CSV = @(data,colnames,rownames,filename) writetable(array2table(data,'VariableNames',colnames,...
                                                                                'RowNames',rownames),...
                                                               filename,'WriteRowNames',true);
writeCell2CSV = @(data,colnames,rownames,filename) writetable(cell2table(data,'VariableNames',colnames,...
                                                                              'RowNames',rownames),...
                                                              filename,'WriteRowNames',true);
makeBndLabels = @(data) arrayfun(@(x) strcat('bnd_pnt_',num2str(x)),1:size(data,1),'UniformOutput',false);
colorscale  = @(data) round(255.*bsxfun(@rdivide, bsxfun(@minus, data, min(data,[],1)), range(data)));
colorscaleg = @(data) round(255.*bsxfun(@rdivide, bsxfun(@minus, data, min(data(:))), range(data(:))));

assignin('caller','writeArray2CSV',writeArray2CSV);
assignin('caller','writeCell2CSV',writeCell2CSV);
assignin('caller','makeBndLabels',makeBndLabels);
assignin('caller','colorscale',colorscale);
assignin('caller','colorscaleg',colorscaleg);
assignin('caller','drawSources',@drawSources);
assignin('caller','drawScatter',@drawScatter);
assignin('caller','drawPortfolioSelections',@drawPortfolioSelections);
assignin('caller','drawPortfolioFootprint',@drawPortfolioFootprint);
assignin('caller','drawGoodBadFootprint',@drawGoodBadFootprint);
assignin('caller','drawFootprint',@drawFootprint);
assignin('caller','drawBinaryPerformance',@drawBinaryPerformance);
assignin('caller','drawBoundary',@drawBoundary);
assignin('caller','getPolygonPoints',@getPolygonPoints);
assignin('caller','footprintBoundary',@footprintBoundary);
assignin('caller','traceAlphaBoundary',@traceAlphaBoundary);
assignin('caller','resolveViewAngle',@resolveViewAngle);

end
% =========================================================================
% SUBFUNCTIONS
% =========================================================================
function lims = axisLimits(Z)
% Interleaved [xmin xmax ymin ymax (zmin zmax)] for axis(...), sized to
% match Z's actual dimensionality (2D or 3D) instead of assuming 2D.
% 'omitnan' guards against instances with an unprojected/missing
% coordinate; without it, a single NaN in Z makes min/max return NaN and
% axis(...) errors instead of just excluding that instance from the bounds.
%
% Padding is 10% of each dimension's own data range rather than a fixed
% +/-1: PILOT's PLS method routinely produces Z values an order of
% magnitude smaller than the standard/analytic methods (e.g. a
% projection spanning +/-0.05), and a fixed +/-1 margin -- plus the
% ceil/floor rounding this used to apply to the bounds themselves --
% swamped that real range entirely, leaving the whole point cloud
% squeezed into a tiny fraction of a much larger empty axes box.
ubound = max(Z,[],1,'omitnan');
lbound = min(Z,[],1,'omitnan');
bad = isnan(ubound) | isnan(lbound);
ubound(bad) = 1;
lbound(bad) = -1;
range_ = ubound - lbound;
range_(range_ == 0) = 1; % all instances share this coordinate; avoid zero padding
pad = 0.1 * range_;
lims = reshape([lbound-pad; ubound+pad], 1, []);
end
% =========================================================================
function labelAxes(is3D)
% Common z1/z2(/z3) axis labelling, so each drawing function doesn't repeat
% the is3D branch just to add zlabel.
xlabel('z_{1}'); ylabel('z_{2}');
if is3D; zlabel('z_{3}'); end
end
% =========================================================================
function plotLine(Z, idx, varargin)
% line(...) dispatches to a 3D plot when given a Z argument; branch once
% here instead of in every caller.
if size(Z,2) == 3
    line(Z(idx,1), Z(idx,2), Z(idx,3), varargin{:});
else
    line(Z(idx,1), Z(idx,2), varargin{:});
end
end
% =========================================================================
function applyView(is3D, viewAngle)
% Applies the optimised camera viewpoint (spec §5.2) if one was supplied,
% otherwise falls back to MATLAB's default isometric 3D view. No-op in 2D.
if ~is3D, return; end
if isempty(viewAngle)
    view(3);
else
    view(viewAngle(1), viewAngle(2));
end
drawCompass(is3D);
end
% =========================================================================
function drawCompass(is3D)
% Coloured z1/z2/z3 axis-direction overlay (spec §8; Figure 5 of Simpson
% et al., 2025), called from applyView so every 3D plot gets it uniformly.
% Anchored at a corner of the current axis limits (already set by the
% caller via axis(axisLimits(Z)) before applyView runs) and drawn in the
% same 3D data coordinates as the scatter/footprint, so it rotates
% together with the rest of the plot instead of staying screen-fixed.
if ~is3D, return; end
lims = axis;
anchor  = [lims(1), lims(3), lims(5)];
armLen  = 0.15 * max([lims(2)-lims(1), lims(4)-lims(3), lims(6)-lims(5)]);
axColors = [0.85 0.10 0.10; 0.10 0.60 0.10; 0.10 0.40 0.85]; % z1/z2/z3
axLabels = {'z_{1}', 'z_{2}', 'z_{3}'};
dirs = eye(3);
% Preserve/restore the caller's hold state rather than forcing it off:
% these helpers are injected into other workspaces via scriptfcn and may
% be composed with further overlays that rely on hold already being on.
ax = gca;
wasHeld = ishold(ax);
hold(ax, 'on');
for i = 1:3
    tip = anchor + armLen * dirs(i,:);
    % HandleVisibility off + explicit Parent: these are a fixed decorative
    % overlay, not plotted data -- without this, a legend() created later
    % in the same axes (AutoUpdate is on by default) would pick up these
    % lines/labels as spurious entries, and they'd be selectable targets
    % for interactive clicks like any other plotted object.
    line([anchor(1) tip(1)], [anchor(2) tip(2)], [anchor(3) tip(3)], ...
        'Color', axColors(i,:), 'LineWidth', 2, ...
        'Parent', ax, 'HandleVisibility', 'off');
    text(tip(1), tip(2), tip(3), axLabels{i}, 'Color', axColors(i,:), ...
        'FontWeight', 'bold', 'FontSize', 10, ...
        'Parent', ax, 'HandleVisibility', 'off');
end
if ~wasHeld
    hold(ax, 'off');
end
end
% =========================================================================
function viewAngle = resolveViewAngle(viewpoint, algoIdx)
% Resolves the [azimuth elevation] (degrees, for view(az,el)) to use for a
% plot, given a model.pilot.viewpoint struct (PILOTviewpoint's output) and
% the algorithm column index the plot is about (or [] for plots that
% aren't tied to a single algorithm, e.g. feature distributions).
%
% Returns [] (meaning: fall back to view(3)) when no viewpoint is
% available, i.e. a 2D projection where PILOTviewpoint was never called.
%
% algoIdx is matched against viewpoint.groups (opts.pilot.viewGroups, or
% the single default group spanning every algorithm); an algorithm not
% covered by any group -- or a plot with no specific algorithm -- falls
% back to the first computed viewpoint.
if isempty(viewpoint)
    viewAngle = [];
    return;
end
g = [];
if ~isempty(algoIdx)
    g = find(cellfun(@(grp) any(grp == algoIdx), viewpoint.groups), 1);
end
if isempty(g)
    g = 1;
end
viewAngle = rad2deg([viewpoint.azimuth(g), viewpoint.elevation(g)]);
end
% =========================================================================
function d = dotDiameter()
% Shared marker diameter (points) for every instance dot drawn by the
% functions below. A single source of truth so drawSources/
% drawPortfolioSelections/drawPortfolioFootprint/drawGoodBadFootprint/
% drawBinaryPerformance (all plotLine()-based, Marker 'o') and
% drawScatter (scatter()-based, sized via dotArea() below) render dots of
% the same apparent size -- tune this one value if the rendered PNGs need
% bigger or smaller dots.
d = 6;
end
% =========================================================================
function a = dotArea()
% scatter()/scatter3()'s size argument is marker area in points^2, not
% diameter, so drawScatter can't just reuse dotDiameter() directly --
% convert to the equal-area circle so its dots match the 'o' markers'
% diameter everywhere else in this file.
d = dotDiameter();
a = pi*(d/2)^2;
end
% =========================================================================
function handle = drawSources(Z, S, viewAngle)
if nargin < 3, viewAngle = []; end
is3D = size(Z,2) == 3;
sourcelabels = cellstr(unique(S));
nsources = length(sourcelabels);
clrs = flipud(lines(nsources));
handle = zeros(nsources,1);
for i=nsources:-1:1
    plotLine(Z, S==sourcelabels{i}, 'LineStyle', 'none', ...
             'Marker', 'o', ...
             'Color', clrs(i,:), ...
             'MarkerFaceColor', clrs(i,:), ...
             'MarkerEdgeColor', 'none', ...
             'MarkerSize', dotDiameter());
    handle(i) = patch([0 0],[0 0], clrs(i,:), 'EdgeColor','none');
end
labelAxes(is3D); title('Sources');
legend(handle, sourcelabels, 'Location', 'NorthEastOutside');
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square; axis(axisLimits(Z)); grid on;
applyView(is3D, viewAngle);

end
% =========================================================================
function handle = drawScatter(Z, X, titlelabel, viewAngle)
if nargin < 4, viewAngle = []; end
is3D = size(Z,2) == 3;
if is3D
    handle = scatter3(Z(:,1), Z(:,2), Z(:,3), dotArea(), X, 'filled');
else
    handle = scatter(Z(:,1), Z(:,2), dotArea(), X, 'filled');
end
caxis([0,1])
labelAxes(is3D); title(titlelabel);
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square; axis(axisLimits(Z)); grid on;
applyView(is3D, viewAngle);
colorbar('EastOutside');

end
% =========================================================================
function drawPortfolioSelections(Z, P, algolabels, titlelabel, viewAngle)
if nargin < 5, viewAngle = []; end
is3D = size(Z,2) == 3;
nalgos = length(algolabels);
algolbls = cell(1,nalgos+1);
h = zeros(1,nalgos+1);
isworthy = sum(bsxfun(@eq, P, 0:nalgos))~=0;
clr = flipud(lines(nalgos+1));
for i=0:nalgos
    if ~isworthy(i+1)
        continue;
    end
    plotLine(Z, P==i, 'LineStyle', 'none', ...
                       'Marker', 'o', ...
                       'Color', clr(i+1,:), ...
                       'MarkerFaceColor', clr(i+1,:), ...
                       'MarkerEdgeColor', 'none', ...
                       'MarkerSize', dotDiameter());
    h(i+1) = patch([0 0],[0 0], clr(i+1,:), 'EdgeColor','none');
    if i==0
        algolbls{i+1} = 'None';
    else
        algolbls{i+1} = strrep(algolabels{i},'_',' ');
    end
end
labelAxes(is3D); title(titlelabel);
legend(h(isworthy), algolbls(isworthy), 'Location', 'NorthEastOutside');
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square; axis(axisLimits(Z)); grid on;
applyView(is3D, viewAngle);

end
% =========================================================================
function h = drawPortfolioFootprint(Z, best, P, algolabels, viewAngle)
if nargin < 5, viewAngle = []; end
is3D = size(Z,2) == 3;
nalgos = length(algolabels);
algolbls = cell(1,nalgos+1);
isworthy = sum(bsxfun(@eq, P, 0:nalgos))~=0;
clr = flipud(lines(nalgos+1));
h = zeros(1,nalgos+1);
for i=0:nalgos
    if ~isworthy(i+1)
        continue;
    end
    plotLine(Z, P==i, 'LineStyle', 'none', ...
                       'Marker', 'o', ...
                       'Color', clr(i+1,:), ...
                       'MarkerFaceColor', clr(i+1,:), ...
                       'MarkerEdgeColor', 'none', ...
                       'MarkerSize', dotDiameter());
    h(i+1) = patch([0 0],[0 0], clr(i+1,:), 'EdgeColor','none');
    if i==0
        algolbls{i+1} = 'None';
    else
        drawFootprint(best{i}, clr(i+1,:), 0.3);
        algolbls{i+1} = strrep(algolabels{i},'_',' ');
    end
end
labelAxes(is3D); title('Portfolio footprints');
legend(h(isworthy), algolbls(isworthy), 'Location', 'NorthEastOutside');
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square; axis(axisLimits(Z)); grid on;
applyView(is3D, viewAngle);

end
% =========================================================================
function h = drawGoodBadFootprint(Z, good, Ybin, titlelabel, viewAngle)
if nargin < 5, viewAngle = []; end
is3D = size(Z,2) == 3;
orange = [1.0 0.6471 0.0];
blue = [0.0 0.0 1.0];
lbls = {'GOOD','BAD'};
h = zeros(1,2);
if any(~Ybin)
    plotLine(Z, ~Ybin, 'LineStyle', 'none', ...
                        'Marker', 'o', ...
                        'Color', orange, ...
                        'MarkerFaceColor', orange, ...
                        'MarkerEdgeColor', 'none', ...
                        'MarkerSize', dotDiameter());
    h(2) = patch([0 0],[0 0], orange, 'EdgeColor','none');
end
if any(Ybin)
    plotLine(Z, Ybin, 'LineStyle', 'none', ...
                       'Marker', 'o', ...
                       'Color', blue, ...
                       'MarkerFaceColor', blue, ...
                       'MarkerEdgeColor', 'none', ...
                       'MarkerSize', dotDiameter());
    h(1) = patch([0 0],[0 0], blue, 'EdgeColor','none');
    drawFootprint(good, blue, 0.3);
end
labelAxes(is3D); title([titlelabel ' Footprints']);
legend(h(h~=0), lbls(h~=0), 'Location', 'NorthEastOutside');
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square; axis(axisLimits(Z)); grid on;
applyView(is3D, viewAngle);

end
% =========================================================================
function handle = drawFootprint(footprint, color, alpha)
% plot() dispatches on footprint.polygon's class: polyshape (2D) or
% alphaShape (3D) both accept FaceColor/EdgeColor/FaceAlpha, so no
% dimension branch is needed here.
hold on;
if isempty(footprint) || isempty(footprint.polygon)
    handle = patch([0 0],[0 0], color, 'EdgeColor','none');
    return
end

handle = plot(footprint.polygon,'FaceColor', color, 'EdgeColor','none', 'FaceAlpha', alpha);
hold off;

end
% =========================================================================
function h = drawBinaryPerformance(Z, Ybin, titlelabel, viewAngle)
if nargin < 4, viewAngle = []; end
is3D = size(Z,2) == 3;
orange = [1.0 0.6471 0.0];
blue = [0.0 0.0 1.0];
lbls = {'GOOD','BAD'};
h = zeros(1,2);
if any(~Ybin)
    h(2) = patch([0 0],[0 0], orange, 'EdgeColor','none');
    plotLine(Z, ~Ybin, 'LineStyle', 'none', ...
                        'Marker', 'o', ...
                        'Color', orange, ...
                        'MarkerFaceColor', orange, ...
                        'MarkerEdgeColor', 'none', ...
                        'MarkerSize', dotDiameter());
end
if any(Ybin)
    h(1) = patch([0 0],[0 0], blue, 'EdgeColor','none');
    plotLine(Z, Ybin, 'LineStyle', 'none', ...
                       'Marker', 'o', ...
                       'Color', blue, ...
                       'MarkerFaceColor', blue, ...
                       'MarkerEdgeColor', 'none', ...
                       'MarkerSize', dotDiameter());
end
labelAxes(is3D); title(titlelabel);
legend(h(h~=0), lbls(h~=0), 'Location', 'NorthEastOutside');
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square; axis(axisLimits(Z)); grid on;
applyView(is3D, viewAngle);

end
% =========================================================================
function handle = drawBoundary(Z, Zedge, titlelabel)
% Draws CLOISTER's empirical space boundary (Zedge, a closed polygon --
% see core/CLOISTER.m) as a red outline over a grey scatter of every
% instance, so the boundary can be read against the point cloud it
% bounds (spec deferred item, #32).
%
% 2D only: CLOISTER's Zedge/Zecorr are computed via a 2D-only convex hull
% (core/CLOISTER.m's convhull(Z(:,1),Z(:,2))) even when the projection
% itself is 3D, so an accurate 3D boundary isn't available yet -- callers
% must not invoke this for a 3D projection (opts.pilot.dims==3); see #32.
handle = scatter(Z(:,1), Z(:,2), dotArea(), [0.7 0.7 0.7], 'filled');
hold on;
idxClosed = [1:size(Zedge,1), 1]; % close the polygon back to its first vertex
plot(Zedge(idxClosed,1), Zedge(idxClosed,2), 'r-', 'LineWidth', 1.5, ...
    'DisplayName', 'CLOISTER empirical bound');
hold off;
labelAxes(false); title(titlelabel);
legend('Location', 'NorthEastOutside');
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square; axis(axisLimits(Z)); grid on;
end
% =========================================================================
function pts = getPolygonPoints(polygon)
% Extract vertex/point matrix from either a polyshape or alphaShape object.
if isa(polygon, 'alphaShape')
    pts = polygon.Points;
else
    pts = polygon.Vertices;
    pts = pts(~any(isnan(pts),2),:); % polyshape uses NaN rows as region separators
end
end
% =========================================================================
function verts = footprintBoundary(fp)
% Extract ordered 2-D boundary vertices from a footprint struct.
% Handles both polyshape (legacy) and alphaShape (TRACE3) polygons.
verts = [];
if ~isfield(fp, 'polygon') || isempty(fp.polygon)
    return;
end
poly = fp.polygon;
if isa(poly, 'alphaShape')
    if size(poly.Points, 2) == 3
        % 3D boundary CSV export not supported; skip silently.
        % Full 3D boundary extraction is deferred to a later phase.
        return;
    end
    verts = traceAlphaBoundary(poly);
elseif isa(poly, 'polyshape')
    verts = poly.Vertices;
end
end
% =========================================================================
function verts = traceAlphaBoundary(poly)
% Trace ordered closed-polygon vertices from a 2-D alphaShape, covering
% EVERY disconnected region (#31) instead of silently returning only the
% first one traced -- the previous version built one adjacency graph from
% boundaryFacets(poly)'s combined edge list (every region's edges
% flattened together) and stopped as soon as it ran out of same-loop
% neighbours, silently dropping every other region. Regions are traced
% independently via boundaryFacets(poly, regionID) instead.
%
% Multiple regions are separated by a row of NaNs in the returned matrix
% -- the same multi-part-polyline convention polyshape.Vertices already
% uses (see getPolygonPoints above): plot()/line() naturally break the
% line at a NaN row, and a consumer that doesn't expect multiple regions
% can split on isnan(verts(:,1)). A single-region shape (the common case)
% is returned identically to before: no separator rows at all.
verts = [];
for r = 1:numRegions(poly)
    [bf, bv] = boundaryFacets(poly, r);
    if isempty(bf)
        continue;
    end
    regionVerts = traceOneRegion(bf, bv);
    if isempty(regionVerts)
        continue;
    end
    if ~isempty(verts)
        verts = [verts; NaN(1, size(regionVerts,2))]; %#ok<AGROW>
    end
    verts = [verts; regionVerts]; %#ok<AGROW>
end
end
% =========================================================================
function verts = traceOneRegion(bf, bv)
% Trace an ordered closed polygon from a single region's boundary-facets
% edge list. bf: (m x 2) edge index pairs into bv's rows (bv may be the
% full alphaShape's point list rather than one scoped to this region
% alone; only the vertices bf actually references matter here). Works
% correctly for a simple, single connected boundary -- traceAlphaBoundary
% calls this once per region rather than once for a whole (possibly
% multi-region) shape, which is what makes that assumption safe again.
regionVertIdx = unique(bf(:));
nRegionVerts = numel(regionVertIdx);
if nRegionVerts == 0, verts = []; return; end
% Build adjacency: each vertex has exactly 2 neighbours on this region's
% boundary. Sized to the largest referenced index (not size(bv,1)) since
% bf's indices may only be a subset of bv's full row range.
maxIdx = max(regionVertIdx);
adj = zeros(maxIdx, 2);
cnt = zeros(maxIdx, 1);
for k = 1:size(bf, 1)
    v1 = bf(k,1); v2 = bf(k,2);
    cnt(v1) = cnt(v1)+1;
    if cnt(v1) <= 2, adj(v1, cnt(v1)) = v2; end
    cnt(v2) = cnt(v2)+1;
    if cnt(v2) <= 2, adj(v2, cnt(v2)) = v1; end
end
% Trace starting from a vertex actually in THIS region (bf(1,1)), not a
% hardcoded global index 1 -- vertex 1 of the whole shape may belong to a
% different region entirely once bf is scoped to region r.
order = zeros(nRegionVerts, 1);
order(1) = bf(1,1);
prev = 0; curr = order(1);
for k = 2:nRegionVerts
    nxt = adj(curr, adj(curr,:) ~= prev & adj(curr,:) ~= 0);
    if isempty(nxt), break; end
    order(k) = nxt(1);
    prev = curr;
    curr = order(k);
end
valid = order ~= 0;
verts = bv(order(valid), :);
end