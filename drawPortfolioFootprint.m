function h = drawPortfolioFootprint(Z, best, Y, algolabels)

% -------------------------------------------------------------------------
% Color definitions
nalgos = length(algolabels);
colors = rand(nalgos,3);
h = zeros(1,nalgos);
for i=1:nalgos
    drawFootprint(best{i}, colors(i,:), 0.2);
    h(i) = line(Z(Y==i,1), Z(Y==i,2), 'LineStyle', 'none', ...
                                      'Marker', '.', ...
                                      'Color', colors(i,:), ...
                                      'MarkerFaceColor', colors(i,:), ...
                                      'MarkerSize', 6);
end
xlabel('z_{1}'); ylabel('z_{2}'); title('Portfolio Footprints');
legend(h, algolabels, 'Location', 'NorthEastOutside');
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square;
end