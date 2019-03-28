function footprint = calculateFootprintPerformance(footprint,Z,Y)
% 
if isempty(footprint.polygon)
    footprint.elements = 0;
    footprint.goodElements = 0;
    footprint.area = 0;
    footprint.density = 0;
    footprint.purity = 0;
    return
end

elements     = false(size(Z,1),1);
goodElements = false(size(Z,1),1);
for i=1:size(footprint.polygon,3)
    elements     = elements     |  inpolygon(Z(:,1),Z(:,2),footprint.polygon(:,1,i),footprint.polygon(:,2,i));
    goodElements = goodElements | (inpolygon(Z(:,1),Z(:,2),footprint.polygon(:,1,i),footprint.polygon(:,2,i)) & ~Y);
end
footprint.elements     = elements;
footprint.goodElements = goodElements;
footprint.area         = sum(footprint.polyArea);
footprint.density      = sum(elements)./footprint.area;
footprint.purity       = sum(goodElements)./sum(elements);
end