function footprint = calculateFootprintPerformance(footprint,x,y)
% 
elements     = false(size(x,1),1);
goodElements = false(size(x,1),1);
for i=1:size(footprint.polygon,3)
    elements     = elements     |  inpolygon(x(:,1),x(:,2),footprint.polygon(:,1,i),footprint.polygon(:,2,i));
    goodElements = goodElements | (inpolygon(x(:,1),x(:,2),footprint.polygon(:,1,i),footprint.polygon(:,2,i)) & ~y);
end
footprint.elements     = elements;
footprint.goodElements = goodElements;
footprint.area         = sum(footprint.polyArea);
footprint.density      = sum(elements)./footprint.area;
footprint.purity       = sum(goodElements)./sum(elements);
end