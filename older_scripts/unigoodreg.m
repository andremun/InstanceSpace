% probgood = mean(model.data.Ybin,1)';
% [pgoodsort,idx] = sort(probgood,'descend');
% model.data.algolabels(idx(pgoodsort>0.1))';
% 
% puniquegood = mean(model.data.Ybin & (sum(model.data.Ybin,2)<=1))';
% [puniquegoodsort,idx] = sort(puniquegood,'descend');
% model.data.algolabels(idx(puniquegoodsort>0.01))';
% 
% 
% numgoodalgos = sum(model.data.Ybin,2);
% percgoodalgos = mean(bsxfun(@eq,numgoodalgos,1:nalgos));
% 
% 
% 
% 
% model = load('C:\Users\mariom1\OneDrive - The University of Melbourne\Documents\MATLAB\InstanceSpace_Regression\trial\results_r1_main_paper_results\model.mat');
% scriptfcn;
clf;
h = drawSources(model.pilot.Z,model.data.S);
line([0 0],[-100 100],'LineStyle', '--', 'color', [0.6 0.6 0.6]);
line([-100 100],[0 0],'LineStyle', '--', 'color', [0.6 0.6 0.6]);
text(3.5, 3,'Q1');
text(-4, 3,'Q2');
text(3.5,-3,'Q4');
text(-4,-3,'Q3');
legend(h, {'BlackBoxData','EvolvedBlackBox','M3C','Repositories'})
% print(gcf,'-dpng',[rootdir 'distribution_sources_quadrants.png']);
% 
% 

rootdir = 'C:\Users\mariom1\OneDrive - The University of Melbourne\Documents\MATLAB\InstanceSpace_Regression\trial\results_r1_main_paper_results\';
h = zeros(1,3);
for ii=1:length(model.data.algolabels)
    clf;
    alpha = 0.4;
    blue = [0.0 0.0 1.0];
    h(1:2) = drawBinaryPerformance(model.pilot.Z, model.data.Ybin(:,ii), ...
                              strrep(model.data.algolabels{ii},'_',' '));
    h(3) = drawFootprint(model.trace.good{ii}, blue, alpha);
    legend(h,{'GOOD','BAD','FTPRN'});
    print(gcf,'-dpng',[rootdir 'footprint_' model.data.algolabels{ii} '.png']);
end


ii = 3; clf; drawBinaryPerformance(model.pilot.Z, presult.Yhat(:,ii), strrep(model.data.algolabels{ii},'_',' '));