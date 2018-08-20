load('E:\InstanceSpace_WorkforceScheduling\rawdata.mat')

opts.thresholds.DIVTHRESHOLD = 0.01; % Minimum percentage allowed of repeated values [?]
opts.thresholds.CORTHRESHOLD = 3;    % Top N features (by correlation) per algorithm that are selected

opts.thresholds.BETATHRESHOLD = 0.5; % Beta-easy threshold

opts.clust.KDEFAULT = 10;            % Default maximum number of clusters
opts.clust.SILTHRESHOLD = 0.55;      % Minimum accepted value for the average silhoute value
opts.clust.NTREES = 50;              % Number of trees for the Random Forest (to determine highest separability in the 2-d projection)
opts.clust.MaxIter = 1000;
opts.clust.Replicates = 100;
opts.clust.UseParallel = false;

opts.footprint.RHO = 10;             % Density threshold
opts.footprint.PI = 0.75;            % Purity threshold
opts.footprint.LOWER_PCTILE = 1;     % 
opts.footprint.UPPER_PCTILE = 25;

opts.pbldr.ntries = 10;
opts.pbldr.analytic = false;
opts.pbldr.cmaopts = bipopcmaes;
opts.pbldr.cmaopts.StopFitness = 0;
opts.pbldr.cmaopts.MaxRestartFunEvals = 0;
opts.pbldr.cmaopts.MaxFunEvals  = 1e4;
opts.pbldr.cmaopts.EvalParallel = 'no';
opts.pbldr.cmaopts.DispFinal = 'off';

out = matilda(X, Y, Ybin, opts);