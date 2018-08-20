function baseFootprint = calculateFootprintCollisions(baseFootprint,testFootprint)
% 
numBasePolygons = size(baseFootprint.polygon,3);
numTestPolygons = size(testFootprint.polygon,3);
totalTestArea   = zeros(numBasePolygons,numTestPolygons);
testArea        = testFootprint.polyArea;
for i=1:numBasePolygons
    basePolygon  = baseFootprint.polygon(:,:,i);
    testPolygon1 = squeeze(testFootprint.polygon(:,1,:));
    testPolygon2 = squeeze(testFootprint.polygon(:,2,:));
    for j=1:numTestPolygons
        xx = polybool('intersection',basePolygon(:,1),basePolygon(:,2),...
                                     testPolygon1(:,j),testPolygon2(:,j));
        if ~isempty(xx), totalTestArea(i,j) = testArea(j); end
    end
end
totalTestArea = sum(totalTestArea,2);
idx = baseFootprint.polyArea<totalTestArea;

baseFootprint.polyArea = baseFootprint.polyArea(~idx);
baseFootprint.polyDensity = baseFootprint.polyDensity(~idx);
baseFootprint.polyElements = baseFootprint.polyElements(~idx);
baseFootprint.polyGoodElements = baseFootprint.polyGoodElements(~idx);
baseFootprint.polyPurity = baseFootprint.polyPurity(~idx);
baseFootprint.polygon = baseFootprint.polygon(:,:,~idx);

end