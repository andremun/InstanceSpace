function model = TRACEtest(model, Z, Ybin, P, beta, algolabels)

nalgos = length(model.best);
disp('-------------------------------------------------------------------------');
disp('  -> TRACE is calculating the algorithm footprints.');
model.test.best = zeros(nalgos,5);
model.test.good = zeros(nalgos,5);
% model.test.bad = zeros(nalgos,5);
% Use the actual data to calculate the footprints
for i=1:nalgos
    model.test.best(i,:) = TRACEtestsummary(model.best{i}, Z,  P==i, model.space.area, model.space.density);
    model.test.good(i,:) = TRACEtestsummary(model.good{i}, Z,  Ybin(:,i), model.space.area, model.space.density);
    % model.test.bad(i,:)  = TRACEtestsummary(model.bad{i},  Z, ~Ybin(:,i), model.space.area, model.space.density);
end

% -------------------------------------------------------------------------
% Beta hard footprints. First step is to calculate them.
disp('-------------------------------------------------------------------------');
disp('  -> TRACE is calculating the beta-footprints.');
% model.test.easy = TRACEtestsummary(model.easy, Z,  beta, model.space.area, model.space.density);
model.test.hard = TRACEtestsummary(model.hard, Z, ~beta, model.space.area, model.space.density);
% -------------------------------------------------------------------------
% Calculating performance
disp('-------------------------------------------------------------------------');
disp('  -> TRACE is preparing the summary table.');
model.summary = cell(nalgos+1,11);
model.summary(1,2:end) = {'Area_Good',...
                          'Area_Good_Normalized',...
                          'Density_Good',...
                          'Density_Good_Normalized',...
                          'Purity_Good',...
                          'Area_Best',...
                          'Area_Best_Normalized',...
                          'Density_Best',...
                          'Density_Best_Normalized',...
                          'Purity_Best'};
model.summary(2:end,1) = algolabels(1:nalgos);
model.summary(2:end,2:end) = num2cell([model.test.good model.test.best]);

disp('  -> TRACE has completed. Footprint analysis results:');
disp(' ');
disp(model.summary);

end
% =========================================================================
function out = TRACEtestsummary(footprint, Z, Ybin, spaceArea, spaceDensity)
% 
if isempty(footprint.polygon) || all(~Ybin)
    out = zeros(5,1);
else
    elements = sum(isinterior(footprint.polygon, Z));
    goodElements = sum(isinterior(footprint.polygon, Z(Ybin,:)));
    density = elements./footprint.area;
    purity = goodElements./elements;
    
    out = [footprint.area,...
           footprint.area/spaceArea,...
           density,...
           density/spaceDensity,...
           purity];
end
out(isnan(out)) = 0;
end
% =========================================================================