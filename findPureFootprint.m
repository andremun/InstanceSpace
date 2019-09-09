function footprint = findPureFootprint(Z, Y, opts)
% If there is no Y to work with, then there is not point on this one
try
    Ig = unique(Z(Y,:),'rows');   % There might be points overlapped, so eliminate them to avoid problems
    numInst = size(Ig,1);
    % We use a random subset of instances to estimate the footprints. Once
    % we do that, we calculate a delaunay triangulation.
    disp(['      -> Using ' num2str(round(100.*opts.PCTILE,1)) ...
          '% of the instances to calculate the footprint.']);
    elim = true(numInst,1);
    aux = cvpartition(numInst,'HoldOut',opts.PCTILE);
    elim(convhull(Ig(:,1),Ig(:,2))) = false; % Keep the ones in the convex hull
    elim(aux.test) = false;
    Ig = Ig(~elim,:);
    disp(['      -> ' num2str(sum(~elim)) ...
          ' instances are going to be used.']);
    disp('      -> Calculating a Delaunay triangulation.');
    polygon = delaunay(Ig);
catch ME
    disp('      -> There is not enough instances to calculate a footprint.');
    disp('      -> Either the algorithm is very weak or the subset of instances used is too small.');
    disp(['ID: ' ME.identifier]);
    footprint.polyArea = [];
    footprint.polyDensity = [];
    footprint.polyElements = 0;
    footprint.polyGoodElements = 0;
    footprint.polyPurity = [];
    footprint.polygon = [];
    footprint.pieces = 0;
    return;
end
numPolygons = size(polygon,1);
% Calculating cardinality, cardinality of good elements, and area of each
% polygon.
polyElements = zeros(numPolygons,1);
polyGoodElements = zeros(numPolygons,1);
polyArea = zeros(numPolygons,1);
for i=1:numPolygons
    polyElements(i) = sum(inpolygon(Z(:,1),Z(:,2),Ig(polygon(i,:),1),Ig(polygon(i,:),2)));
    polyGoodElements(i) = sum(inpolygon(Z(Y,1),Z(Y,2),Ig(polygon(i,:),1),Ig(polygon(i,:),2)));
    polyArea(i) = polyarea(Ig(polygon(i,:),1),Ig(polygon(i,:),2));
end
polyDensity = polyElements./polyArea;
polyPurity = polyGoodElements./polyElements;
% Which polygons fulfill the density and purity requriements
idx = (polyDensity>=opts.RHO) & (polyPurity>=opts.PI);
% Eliminate any of these poligons. If they don't fulfill the density/purity
% requirements, they only make the footprint worse.
polyArea    = polyArea(idx);
polyDensity = polyDensity(idx);
polyElements = polyElements(idx);
polyGoodElements = polyGoodElements(idx);
polyPurity = polyPurity(idx);
polygon = polygon(idx,:);
numPolygons = size(polygon,1);
disp(['      -> ' num2str(round(100.*mean(idx),1)) '%' ...
      ' of the footprint fullfils the density and purity requirements.']);
if mean(Y)<0.98
    % Now, sort the elements from largest to lowest and let's see if we can
    % improve the density/purity by removing the worst 25%
    [~,idx1] = sort(polyArea,'descend'); % big to small
    [~,idx2] = sort(polyDensity,'ascend'); % bad to worse
    [~,idx3] = sort(polyPurity,'ascend');
    worst = round(0.02.*numPolygons); % the outliers
    aux = idx1(any(bsxfun(@eq,idx1(1:worst),idx2(1:worst)') | ...
                   bsxfun(@eq,idx1(1:worst),idx3(1:worst)'),2));
    idx = true(numPolygons,1);
    idx(aux) = false;
%     footDensity = sum(polyElements)./sum(polyArea);
%     footPurity = sum(polyGoodElements)./sum(polyElements);
%     for i=1:numPolygons
%         % If this polygon does not meet the density/purity requirements should it be discarded?
%         if ~idx(i)
%             sumCardPoly = sum(polyElements(idx));
%             densityWoutPoly = sumCardPoly./sum(polyArea(idx));
%             purityWoutPoly = sum(polyGoodElements(idx))./sumCardPoly;
%             if ((footDensity>densityWoutPoly) && (footPurity>purityWoutPoly)) || ...
%                ((opts.RHO>densityWoutPoly) && (opts.PI>purityWoutPoly))     %#ok<BDSCI>
%                 % If the footprint no longer meets the density/purity
%                 % requirements without the polygon, keep it.
%                 idx(i) = true;
%                 footDensity = sum(polyElements(idx))./sum(polyArea(idx));
%                 footPurity = sum(polyGoodElements(idx))./sum(polyElements(idx));
%             end
%         end
%     end
    polygon     = polygon(idx,:);
    numPolygons = size(polygon,1);
    % Setting the footprint data
    footprint.polyArea    = polyArea(idx);
    footprint.polyDensity = polyDensity(idx);
    footprint.polyElements = polyElements(idx);
    footprint.polyGoodElements = polyGoodElements(idx);
    footprint.polyPurity = polyPurity(idx);
else
    % Setting the footprint data
    footprint.polyArea    = polyArea;
    footprint.polyDensity = polyDensity;
    footprint.polyElements = polyElements;
    footprint.polyGoodElements = polyGoodElements;
    footprint.polyPurity = polyPurity;
end

footprint.polygon = zeros(3,2,numPolygons);
% graph       = zeros(numInst);
for i=1:numPolygons
    [polyx,polyy] = poly2cw(Ig(polygon(i,:),1),Ig(polygon(i,:),2));
    footprint.polygon(:,:,i) = [polyx polyy];
%     graph(polygon(i,:),polygon(i,:)) = squareform(pdist(Ig(polygon(i,:),:)));
end
% graph       = sparse(graph);
% distGraph   = zeros(numInst);
% for i=1:numInst
%     distGraph(i,:) = graphshortestpath(graph,i);
% end
% finish = find(~isinf(distGraph(1,:)),1,'last');
% start  = finish + 1;
% footprint.pieces = 1;
% while start<numInst && finish<=numInst
%     finish = start + find(~isinf(distGraph(start,start:numInst)),1,'last');
%     start  = finish + 1;
%     footprint.pieces = footprint.pieces + 1;
% end

end