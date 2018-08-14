function out = calculateFootprintPerformance(footprint, spaceArea, spaceDensity)
% 
aux = sum(footprint.elements);
out = zeros(1,5);
out(1) = sum(footprint.area);
out(2) = out(1)/spaceArea;
out(3) = aux./out(1);
out(4) = out(3)/spaceDensity;
out(5) = sum(footprint.goodElements)./aux;
end