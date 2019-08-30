function drawSVMPortfolioSelections(Z, psel, algolabels)

nalgos = length(algolabels);
isworty = sum(bsxfun(@eq, psel, 1:nalgos))~=0;
clr = parula(nalgos);
mrkrs = {'o','s','x'};
cnt = 1;
for i=1:nalgos
    if ~isworty(i)
        continue;
    end
    line(Z(psel==i,1), Z(psel==i,2), 'LineStyle', 'none', ...
                                     'Marker', mrkrs{cnt}, ...
                                     'Color', clr(i,:), ...
                                     'MarkerFaceColor', clr(i,:), ...
                                     'MarkerSize', 6);
    cnt = cnt+1;
    if cnt>length(mrkrs)
        cnt = 1;
    end
end
xlabel('z_{1}'); ylabel('z_{2}');
legend(algolabels(isworty), 'Location', 'NorthEastOutside');
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square;

end