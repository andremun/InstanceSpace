function drawSVMPortfolioSelections(Z, psel, algolabels)

nalgos = length(algolabels);
algolbls = cell(1,nalgos+1);
isworty = sum(bsxfun(@eq, psel, 0:nalgos))~=0;
clr = parula(nalgos+1);
mrkrs = {'.','+','x'};
cnt = 1;
for i=0:nalgos
    if ~isworty(i+1)
        continue;
    end
    line(Z(psel==i,1), Z(psel==i,2), 'LineStyle', 'none', ...
                                     'Marker', mrkrs{cnt}, ...
                                     'Color', clr(i+1,:), ...
                                     'MarkerFaceColor', clr(i+1,:), ...
                                     'MarkerSize', 6);
    cnt = cnt+1;
    if cnt>length(mrkrs)
        cnt = 1;
    end
    if i==0
        algolbls{i+1} = 'None';
    else
        algolbls{i+1} = strrep(algolabels{i},'_',' ');
    end
end
xlabel('z_{1}'); ylabel('z_{2}');
legend(algolbls(isworty), 'Location', 'NorthEastOutside');
set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(findall(gcf,'-property','LineWidth'),'LineWidth',1);
axis square;

end