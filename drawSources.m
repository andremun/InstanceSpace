function h = drawSources(Z, S)

sourcelabels = cellstr(unique(S));
nsources = length(sourcelabels);
clrs = flipud(jet(nsources+1));
h = zeros(nsources,1);
for i=nsources:-1:1
    h(i) = line(Z(S==sourcelabels{i},1), ...
                Z(S==sourcelabels{i},2), ...
                'LineStyle', 'none', ...
                'Marker', '.', ...
                'Color', clrs(i,:), ...
                'MarkerFaceColor', clrs(i,:), ...
                'MarkerSize', 6);
end
xlabel('z_{1}'); ylabel('z_{2}'); title('Sources');
legend(h, sourcelabels, 'Location', 'NorthEastOutside');
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square;

end