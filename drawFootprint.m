function handle = drawFootprint(footprint, color)
% 
hold on;
for i=1:size(footprint.polygon,3)
    aux = patch(footprint.polygon(:,1,i), ...
                footprint.polygon(:,2,i), ...
                color, ...
                'EdgeColor','none');
    if i==1, handle = aux; end
end
hold off;
end