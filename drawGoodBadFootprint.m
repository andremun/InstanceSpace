function h = drawGoodBadFootprint(Z, Ybin, good, titlelabel)

clf;
lbls = {'GOOD','BAD'};
h = zeros(1,2);
if any(~Ybin)
    h(2) = line(Z(~Ybin,1), Z(~Ybin,2), 'LineStyle', 'none', ...
                                     'Marker', '.', ...
                                     'Color', [0.8 0.8 0.8], ...
                                     'MarkerFaceColor', [0.8 0.8 0.8], ...
                                     'MarkerSize', 6);
end
if any(Ybin)
    h(1) = drawFootprint(good, [0.0 0.0 0.0]);     % dark green
    line(Z(Ybin,1), Z(Ybin,2), 'LineStyle', 'none', ...
                               'Marker', '.', ...
                               'Color', [0.0 0.0 0.0], ...
                               'MarkerFaceColor', [0.0 0.0 0.0], ...
                               'MarkerSize', 6);
end
% handle(2) = drawFootprint(bad, [1.0 0.0 0.0]);      % Red
xlabel('z_{1}'); ylabel('z_{2}'); title(titlelabel);
legend(h(h~=0), lbls(h~=0), 'Location', 'NorthEastOutside');
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square;
end