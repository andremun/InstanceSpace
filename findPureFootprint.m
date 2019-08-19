function footprint = findPureFootprint(Z, Y, opts)
% If there is no Y to work with, then there is not point on this one
try
    Ig = unique(Z(Y,:),'rows');   % There might be points overlapped, so eliminate them to avoid problems
    numInst = size(Ig,1);
    % ---------
    % We are trying here a different approach
    elim = true(numInst,1);
    aux = cvpartition(numInst,'HoldOut',opts.PCTILE);
    elim(convhull(Ig(:,1),Ig(:,2))) = false; % Keep the ones in the convex hull
    elim(aux.test) = false;
    Ig = Ig(~elim,:);
    polygon = delaunay(Ig);
catch ME
    disp('   -> There is not enough instances to calculate the footprint.');
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
% % ------
% % Eliminate those candidates that are too small or too big
% D         = pdist(Ig)';               % Calculate the distance among all points
% Bounds    = prctile(D,[opts.LOWER_PCTILE opts.UPPER_PCTILE]);
% idx       = D<Bounds(1);
% [row,col] = meshgrid(1:numInst);
% row       = row(tril(true(numInst),-1));
% col       = col(tril(true(numInst),-1));
% idx       = [row(idx) col(idx)];
% elim      = false(numInst,1);
% for i=1:size(idx,1)
%     if ~elim(idx(i,1))
%         elim(idx(i,2)) = true;
%     else
%         elim(idx(i,2)) = false;
%     end
% end
% Ig          = Ig(~elim,:);
% numInst     = size(Ig,1);
% if numInst>3
%     polygon     = delaunay(Ig);
% else
%     footprint.polyArea = [];
%     footprint.polyDensity = [];
%     footprint.polyElements = 0;
%     footprint.polyGoodElements = 0;
%     footprint.polyPurity = [];
%     footprint.polygon = [];
%     footprint.pieces = 0;
%     return;
% end
% numPolygons = size(polygon,1);
% % ------
% elim        = false(numPolygons,1);
% for i=1:numPolygons
%     elim(i) =  any(pdist(Ig(polygon(i,:),:))>Bounds(2));
% end
% polygon     = polygon(~elim,:);
% numPolygons = size(polygon,1);
% % ------
% Calculating cardinality, cardinality of good elements, and area of each
% polygon.
cardPoly = zeros(numPolygons,1);
cardPolyG = zeros(numPolygons,1);
polyArea = zeros(numPolygons,1);
for i=1:numPolygons
    cardPoly(i) = sum(inpolygon(Z(:,1),Z(:,2),Ig(polygon(i,:),1),Ig(polygon(i,:),2)));
    cardPolyG(i) = sum(inpolygon(Z(Y,1),Z(Y,2),Ig(polygon(i,:),1),Ig(polygon(i,:),2)));
    polyArea(i) = polyarea(Ig(polygon(i,:),1),Ig(polygon(i,:),2));
end
polyDensity = cardPoly./polyArea;
polyPurity = cardPolyG./cardPoly;
% Which polygons fulfill the density and purity requriements
idx = (polyDensity>=opts.RHO) & (polyPurity>=opts.PI);
for i=1:numPolygons
    % If this polygon does not meet the density/purity requirements should it be discarded?
    if ~idx(i)
        sumCardPoly = sum(cardPoly(idx));
        densityWoutPoly = sumCardPoly./sum(polyArea(idx));
        purityWoutPoly = sum(cardPolyG(idx))./sumCardPoly;
        if (densityWoutPoly<opts.RHO) && (purityWoutPoly<opts.PI) %#ok<BDSCI>
            % If the footprint no longer meets the density/purity
            % requirements without the polygon, keep it.
            idx(i) = true;
        end
    end
end
polygon     = polygon(idx,:);
numPolygons = size(polygon,1);
% Setting the footprint data
footprint.polyArea    = polyArea(idx);
footprint.polyDensity = polyDensity(idx);
footprint.polyElements = cardPoly(idx);
footprint.polyGoodElements = cardPolyG(idx);
footprint.polyPurity = polyPurity(idx);
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