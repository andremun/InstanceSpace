function baseFootprint = calculateFootprintCollisions(baseFootprint,testFootprint,opts)
% 
if isempty(baseFootprint.polyArea) || isempty(testFootprint.polyArea)
    return;
end

numBasePolygons = size(baseFootprint.polygon,3);
numTestPolygons = size(testFootprint.polygon,3);
idx = true(numBasePolygons,1);
testPolygon1 = squeeze(testFootprint.polygon(:,1,:));
testPolygon2 = squeeze(testFootprint.polygon(:,2,:));
testPolyDensity = testFootprint.polyDensity;
testPolyPurity = testFootprint.polyPurity;

for i=1:numBasePolygons
    basePolygon = polyshape(baseFootprint.polygon(:,1,i),baseFootprint.polygon(:,2,i),'Simplify',false);
    basePolyDensity = baseFootprint.polyDensity(i);
    basePolyPurity = baseFootprint.polyPurity(i);
    for j=1:numTestPolygons
        aux = intersect(basePolygon, polyshape(testPolygon1(:,j),testPolygon2(:,j),'Simplify',false));
        if aux.NumRegions~=0 && (basePolyDensity < testPolyDensity(j) || ...
                                 basePolyPurity < testPolyPurity(j))
            idx(i) = false; % Remove this polygon
        end
    end
end

for i=1:numBasePolygons
    % If this polygon does not meet the density/purity requirements should it be discarded?
    if ~idx(i)
        sumCardPoly = sum(baseFootprint.polyElements(idx));
        densityWoutPoly = sumCardPoly./sum(baseFootprint.polyArea(idx));
        purityWoutPoly = sum(baseFootprint.polyGoodElements(idx))./sumCardPoly;
        if (densityWoutPoly<opts.RHO) && (purityWoutPoly<opts.PI) %#ok<BDSCI>
            % If the footprint no longer meets the density/purity
            % requirements without the polygon, keep it.
            idx(i) = true;
        end
    end
end

baseFootprint.polyArea = baseFootprint.polyArea(idx);
baseFootprint.polyDensity = baseFootprint.polyDensity(idx);
baseFootprint.polyElements = baseFootprint.polyElements(idx);
baseFootprint.polyGoodElements = baseFootprint.polyGoodElements(idx);
baseFootprint.polyPurity = baseFootprint.polyPurity(idx);
baseFootprint.polygon = baseFootprint.polygon(:,:,idx);

end