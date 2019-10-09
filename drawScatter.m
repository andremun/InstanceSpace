function handle = drawScatter(Z, X, titlelabel)

handle = scatter(Z(:,1), Z(:,2), 12, X, 'filled');
caxis([0,1])
xlabel('z_{1}'); ylabel('z_{2}'); title(titlelabel);
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square; colorbar('EastOutside');
end