function out = calculateFootprintSummary(footprint, spaceArea, spaceDensity)
% 
out = zeros(1,5);
out(1) = footprint.area;
out(2) = footprint.area/spaceArea;
out(3) = sum(footprint.elements)./footprint.area;
out(4) = out(3)/spaceDensity;
out(5) = sum(footprint.goodElements)./sum(footprint.elements);
end