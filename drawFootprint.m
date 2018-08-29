function handle = drawFootprint(footprint, color)
% 
hold on;
if isempty(footprint.polygon)
    handle = patch([0 0],[0 0],...
                   color, ...
                   'EdgeColor','none');
    return
end

for i=1:size(footprint.polygon,3)
    aux = patch(footprint.polygon(:,1,i), ...
                footprint.polygon(:,2,i), ...
                color, ...
                'EdgeColor','none');
    if i==1, handle = aux; end
end
hold off;
end