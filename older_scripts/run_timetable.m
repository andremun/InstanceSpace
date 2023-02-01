% ------ Rephrasing the Timetable data
rdir = '../../InstanceSpace_Timetabling/';
load([rdir 'workspace.mat']);
% Loading the options
model.opts.perf.MaxPerf = opts.perf.MaxMin;              % True if Y is a performance measure to maximize, False if it is a cost measure to minimise.
model.opts.perf.AbsPerf = opts.perf.AbsPerf;               % True if an absolute performance measure, False if a relative performance measure
model.opts.perf.epsilon = opts.perf.epsilon;               % Threshold of good performance

model.opts.general = opts.general;      % Beta-easy threshold
model.opts.auto = opts.auto;               % Automatic preprocessing on. Set to false if you don't want any preprocessing
model.opts.bound = opts.bound;                 % Bound the outliers. True if you want to bound the outliers, false if you don't
model.opts.norm = opts.norm;                  % Normalize/Standarize the data. True if you want to apply Box-Cox and Z transformations to stabilize the variance and scale N(0,1)
model.opts.corr = opts.corr;
model.opts.clust = opts.clust;

model.opts.pilot = opts.pbldr;            % Calculate the analytical or numerical solution
model.opts.pythia.cvfolds = opts.oracle.cvfolds;
model.opts.pythia.useweights = false;
model.opts.pythia.uselibsvm = true;
model.opts.pythia.params = out.algosel.paramgrid(out.algosel.paramidx,:);

model.opts.cloister.pval = 0.05;
model.opts.cloister.cthres = 0.7;

model.opts.trace.usesim = false;          % Use the actual or simulated data to calculate the footprints
model.opts.trace.PI = 0.55;               % Purity threshold

model.opts.selvars = opts.selvars;

model.opts.outputs.csv = true;               %
model.opts.outputs.web = true;              % NOTE: MAKE THIS FALSE IF YOU ARE USING THIS CODE LOCALY - This flag is only useful if the system is being used 'online' through matilda.unimelb.edu.au
model.opts.outputs.png = true;   

[model.clust.selvars,idx] = sort(out.clust.selvars);

model.data.X = X(:,idx);
model.data.Y = Y;
model.data.Xraw = Xraw;
model.data.Yraw = Yraw;
model.data.Ybin = Ybin;
model.data.S = S;
model.data.beta = beta;
model.data.bestPerformace = bestPerformace;
model.data.P = portfolio;
model.data.featlabels = featlabels(:,idx);
model.data.algolabels = algolabels;
model.data.instlabels = instlabels;
model.data.numGoodAlgos = sum(model.data.Ybin,2);

model.pilot = out.pbldr;
model.bound = out.bound;
model.norm  = out.norm;
model.corr  = out.corr;
model.clust = out.clust;

% Have to figure out a way to reconstruct the alpha and X0 values
model.pilot.A = model.pilot.A(:,idx);
model.pilot.B = model.pilot.B(idx,:);

aux = [model.pilot.A(:); model.pilot.B(:); model.pilot.C(:)];
alpha = 0.*model.pilot.alpha;
X0 = alpha;
for i=1:length(alpha)
    [row,col] = find(model.pilot.alpha == aux(i));
    alpha(i,:) = model.pilot.alpha(row,:);
    X0(i,:) = model.pilot.X0(row,:);
end
model.pilot.alpha = alpha;
model.pilot.X0 = X0;

auxclust = false(1,7);
auxclust(out.clust.selvars) = true;
aux = false(1,11);
aux(out.corr.selvars) = auxclust;

model.featsel.idx = aux;

save([rdir 'model.mat'],'-struct','model'); % Save the main results

%% 
clearvars;
rdir = '../../InstanceSpace_Timetabling/';
model = load([rdir 'model.mat']);

% opts = model.opts;
% model.pilot.summary = cell(3, length(model.data.featlabels)+1);
% model.pilot.summary(1,2:end) = model.data.featlabels;
% model.pilot.summary(2:end,1) = {'Z_{1}','Z_{2}'};
% model.pilot.summary(2:end,2:end) = num2cell(round(model.pilot.A,4));
% model.cloist = CLOISTER(model.data.X, model.pilot.A, opts.cloister);
% model.pythia = PYTHIA(model.pilot.Z, model.data.Yraw, model.data.Ybin, model.data.bestPerformace, model.data.algolabels, opts.pythia);
% model.trace = TRACE(model.pilot.Z, model.data.Ybin, model.data.P, model.data.beta, model.data.algolabels, opts.trace);
% save([rdir 'model.mat'],'-struct','model'); % Save the main results
% scriptcsv(model,rdir);
scriptpng(model,rdir);
h = drawSources(model.pilot.Z, model.data.S);
nsources = length(h);
h(nsources+1) = line(model.cloist.Zedge(:,1),model.cloist.Zedge(:,2), 'Color', [0.49 0.18 0.56]);
sourcelabels = cellstr(unique(model.data.S));
sourcelabels{nsources+1} = 'estimated bound';
legend(h, sourcelabels, 'Location', 'NorthEastOutside');
axis square; axis([-10 10 -12 12]);
print(gcf,'-dpng',[rdir 'distribution_sources.png']);

%%
clearvars;
rdir = 'C:\Users\mariom1\OneDrive - The University of Melbourne\Documents\MATLAB\InstanceSpace_Timetabling\';
model = load([rdir 'model.mat']);
PI=0:0.05:0.95;
nalgos = length(model.data.algolabels);
summaryExp = zeros(nalgos,10,length(PI));
summarySim = zeros(nalgos,10,length(PI));
for i=1:length(PI)
    opts.trace.PI = PI(i);
    model.trace.exp{i} = TRACE(model.pilot.Z, model.data.Ybin, model.data.P, model.data.beta, model.data.algolabels, opts.trace);
    model.trace.sim{i} = TRACE(model.pilot.Z, model.pythia.Yhat, model.pythia.selection0, model.data.beta, model.data.algolabels, opts.trace);
    summaryExp(:,:,i) = cell2mat(model.trace.exp{i}.summary(2:end,2:end));
    summarySim(:,:,i) = cell2mat(model.trace.sim{i}.summary(2:end,2:end));
end
summary{1} = [PI' squeeze(mean(summaryExp(:,[2 4 5 7 9 10],:),1))'];
summary{2} = [PI' squeeze(mean(summarySim(:,[2 4 5 7 9 10],:),1))'];

