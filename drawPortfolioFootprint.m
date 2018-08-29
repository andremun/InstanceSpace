function handle = drawPortfolioFootprint(best, algolabels)

% -------------------------------------------------------------------------
% Color definitions
nalgos = length(algolabels);
colors = rand(nalgos,3);
handle = zeros(1,nalgos);
for i=1:nalgos
    handle(i) = drawFootprint(best{i}, colors(i,:));
end
xlabel('z_{1}'); ylabel('z_{2}'); title('PORTFOLIO - PERFORMANCE')
legend(handle, algolabels, 'Location', 'NorthEastOutside');
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square;
end