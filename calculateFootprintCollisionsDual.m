function [base,test] = calculateFootprintCollisionsDual(base,test)
% 
if isempty(base.polyArea) || isempty(test.polyArea)
    return;
end

numBasePoly = size(base.polygon,3);
numTestPoly = size(test.polygon,3);
idxBase = true(numBasePoly,1);
idxTest = true(numTestPoly,1);
testPoly1 = squeeze(test.polygon(:,1,:));
testPoly2 = squeeze(test.polygon(:,2,:));
% testPolyArea = test.polyArea;
testPolyDensity = test.polyDensity;
testPolyPurity = test.polyPurity;

for i=1:numBasePoly
    basePoly = polyshape(base.polygon(:,1,i),base.polygon(:,2,i),'Simplify',false);
    % basePolyArea = base.polyArea(i);
    basePolyDensity = base.polyDensity(i);
    basePolyPurity = base.polyPurity(i);
    for j=1:numTestPoly
        if ~idxTest(j)
            continue;
        end
        
        aux = intersect(basePoly, polyshape(testPoly1(:,j),testPoly2(:,j),'Simplify',false));
        if aux.NumRegions~=0 && (basePolyDensity < testPolyDensity(j) || ...
                                 basePolyPurity < testPolyPurity(j))% || ...
                                 %basePolyArea < testPolyArea(j))
            idxBase(i) = false; % Remove this polygon
            break;
        elseif aux.NumRegions~=0 && (basePolyDensity > testPolyDensity(j) || ...
                                     basePolyPurity > testPolyPurity(j))% || ...
                                     %basePolyArea > testPolyArea(j))
            idxTest(j) = false; % Remove this polygon
        end
    end
end

% for i=1:numBasePoly
%     % If this polygon does not meet the density/purity requirements should it be discarded?
%     if ~idxBase(i)
%         sumCardPoly = sum(base.polyElements(idxBase));
%         densityWoutPoly = sumCardPoly./sum(base.polyArea(idxBase));
%         purityWoutPoly = sum(base.polyGoodElements(idxBase))./sumCardPoly;
%         if (densityWoutPoly<opts.RHO) && (purityWoutPoly<opts.PI) %#ok<BDSCI>
%             % If the footprint no longer meets the density/purity
%             % requirements without the polygon, keep it.
%             idxBase(i) = true;
%         end
%     end
% end

disp(['      -> ' num2str(round(100.*mean(idxBase),1)) '%' ...
      ' of the base footprint has no contradictions.']);

base.polyArea = base.polyArea(idxBase);
base.polyDensity = base.polyDensity(idxBase);
base.polyElements = base.polyElements(idxBase);
base.polyGoodElements = base.polyGoodElements(idxBase);
base.polyPurity = base.polyPurity(idxBase);
base.polygon = base.polygon(:,:,idxBase);

% for j=1:numTestPoly
%     % If this polygon does not meet the density/purity requirements should it be discarded?
%     if ~idxTest(j)
%         sumCardPoly = sum(test.polyElements(idxTest));
%         densityWoutPoly = sumCardPoly./sum(test.polyArea(idxTest));
%         purityWoutPoly = sum(test.polyGoodElements(idxTest))./sumCardPoly;
%         if (densityWoutPoly<opts.RHO) && (purityWoutPoly<opts.PI) %#ok<BDSCI>
%             % If the footprint no longer meets the density/purity
%             % requirements without the polygon, keep it.
%             idxTest(j) = true;
%         end
%     end
% end

disp(['      -> ' num2str(round(100.*mean(idxTest),1)) '%' ...
      ' of the test footprint has no contradictions.']);

test.polyArea = test.polyArea(idxTest);
test.polyDensity = test.polyDensity(idxTest);
test.polyElements = test.polyElements(idxTest);
test.polyGoodElements = test.polyGoodElements(idxTest);
test.polyPurity = test.polyPurity(idxTest);
test.polygon = test.polygon(:,:,idxTest);

end