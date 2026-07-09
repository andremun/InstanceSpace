function scriptfcn
% -------------------------------------------------------------------------
% scriptfcn.m
% -------------------------------------------------------------------------
%
% By: Mario Andres Munoz Acosta
%     School of Mathematics and Statistics
%     The University of Melbourne
%     Australia
%     2020
%
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
assignin('caller','getPolygonPoints',@getPolygonPoints);

end
% =========================================================================
% SUBFUNCTIONS
% =========================================================================
function lims = axisLimits(Z)
ubound = ceil(max(Z,[],1,'omitnan'));
lbound = floor(min(Z,[],1,'omitnan'));
lims = reshape([lbound-1; ubound+1], 1, []);
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
function handle = drawSources(Z, S)

is3D = size(Z,2) == 3;
sourcelabels = cellstr(unique(S));
nsources = length(sourcelabels);
clrs = flipud(lines(nsources));
handle = zeros(nsources,1);
for i=nsources:-1:1
    plotLine(Z, S==sourcelabels{i}, 'LineStyle', 'none', ...
             'Marker', '.', ...
             'Color', clrs(i,:), ...
             'MarkerFaceColor', clrs(i,:), ...
             'MarkerSize', 8);
    handle(i) = patch([0 0],[0 0], clrs(i,:), 'EdgeColor','none');
end
labelAxes(is3D); title('Sources');
legend(handle, sourcelabels, 'Location', 'NorthEastOutside');
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square; axis(axisLimits(Z));
if is3D; view(3); end

end
% =========================================================================
function handle = drawScatter(Z, X, titlelabel)

is3D = size(Z,2) == 3;
if is3D
    handle = scatter3(Z(:,1), Z(:,2), Z(:,3), 8, X, 'filled');
else
    handle = scatter(Z(:,1), Z(:,2), 8, X, 'filled');
end
caxis([0,1])
labelAxes(is3D); title(titlelabel);
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square; axis(axisLimits(Z));
if is3D; view(3); end
colorbar('EastOutside');

end
% =========================================================================
function drawPortfolioSelections(Z, P, algolabels, titlelabel)

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
                       'Marker', '.', ...
                       'Color', clr(i+1,:), ...
                       'MarkerFaceColor', clr(i+1,:), ...
                       'MarkerSize', 4);
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
axis square; axis(axisLimits(Z));
if is3D; view(3); end

end
% =========================================================================
function h = drawPortfolioFootprint(Z, best, P, algolabels)

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
                       'Marker', '.', ...
                       'Color', clr(i+1,:), ...
                       'MarkerFaceColor', clr(i+1,:), ...
                       'MarkerSize', 4);
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
axis square; axis(axisLimits(Z));
if is3D; view(3); end

end
% =========================================================================
function h = drawGoodBadFootprint(Z, good, Ybin, titlelabel)

is3D = size(Z,2) == 3;
orange = [1.0 0.6471 0.0];
blue = [0.0 0.0 1.0];
lbls = {'GOOD','BAD'};
h = zeros(1,2);
if any(~Ybin)
    % drawFootprint(bad, orange, 0.2);
    plotLine(Z, ~Ybin, 'LineStyle', 'none', ...
                        'Marker', '.', ...
                        'Color', orange, ...
                        'MarkerFaceColor', orange, ...
                        'MarkerSize', 4);
    h(2) = patch([0 0],[0 0], orange, 'EdgeColor','none');
end
if any(Ybin)
    plotLine(Z, Ybin, 'LineStyle', 'none', ...
                       'Marker', '.', ...
                       'Color', blue, ...
                       'MarkerFaceColor', blue, ...
                       'MarkerSize', 4);
    h(1) = patch([0 0],[0 0], blue, 'EdgeColor','none');
    drawFootprint(good, blue, 0.3);
end
labelAxes(is3D); title([titlelabel ' Footprints']);
legend(h(h~=0), lbls(h~=0), 'Location', 'NorthEastOutside');
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square; axis(axisLimits(Z));
if is3D; view(3); end

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
function h = drawBinaryPerformance(Z, Ybin, titlelabel)

is3D = size(Z,2) == 3;
orange = [1.0 0.6471 0.0];
blue = [0.0 0.0 1.0];
lbls = {'GOOD','BAD'};
h = zeros(1,2);
if any(~Ybin)
    h(2) = patch([0 0],[0 0], orange, 'EdgeColor','none');
    plotLine(Z, ~Ybin, 'LineStyle', 'none', ...
                        'Marker', '.', ...
                        'Color', orange, ...
                        'MarkerFaceColor', orange, ...
                        'MarkerSize', 4);
end
if any(Ybin)
    h(1) = patch([0 0],[0 0], blue, 'EdgeColor','none');
    plotLine(Z, Ybin, 'LineStyle', 'none', ...
                       'Marker', '.', ...
                       'Color', blue, ...
                       'MarkerFaceColor', blue, ...
                       'MarkerSize', 4);
end
labelAxes(is3D); title(titlelabel);
legend(h(h~=0), lbls(h~=0), 'Location', 'NorthEastOutside');
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square; axis(axisLimits(Z));
if is3D; view(3); end

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