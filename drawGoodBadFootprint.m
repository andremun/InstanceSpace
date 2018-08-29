function handle = drawGoodBadFootprint(good, bad, titlelabel)

clf;
handle = zeros(1,2);
handle(1) = drawFootprint(good, [0.0 0.5 0.0]);     % dark green
handle(2) = drawFootprint(bad, [1.0 0.0 0.0]);      % Red
xlabel('z_{1}'); ylabel('z_{2}'); title(titlelabel);
legend(handle,{'GOOD','BAD'}, 'Location', 'NorthEastOutside');
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square;
end