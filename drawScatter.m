function handle = drawScatter(X, Z, titlelabel)

clf;
scale = (X-min(X))./range(X);
handle = scatter(Z(:,1),Z(:,2),14,scale,'filled');
xlabel('z_{1}'); ylabel('z_{2}'); title(titlelabel);
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square; colorbar('EastOutside');
end