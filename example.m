% load('E:\InstanceSpace_WorkforceScheduling\rawdata.mat')
load('E:\InstanceSpace_Knapsack\MATLAB\rawdata.mat');

opts.general.performanceFlag = false;   % True if Y is a performance measure, False if it is a cost measure.
opts.general.betaThreshold = 0.5;       % Beta-easy threshold

opts.diversity.flag = true;             % Run diversity calculation
opts.diversity.threshold = 0.01;        % Minimum percentage allowed of repeated values [?]

opts.corr.flag = true;
% opts.corr.threshold = 3;                % Top N features (by correlation) per algorithm that are selected
opts.corr.threshold = 5;                % Top N features (by correlation) per algorithm that are selected

opts.clust.flag = true;
opts.clust.KDEFAULT = 10;               % Default maximum number of clusters
% opts.clust.SILTHRESHOLD = 0.55;         % Minimum accepted value for the average silhoute value
opts.clust.SILTHRESHOLD = 0.70;         % Minimum accepted value for the average silhoute value
opts.clust.NTREES = 50;                 % Number of trees for the Random Forest (to determine highest separability in the 2-d projection)
opts.clust.MaxIter = 1000;
opts.clust.Replicates = 100;
opts.clust.UseParallel = false;

opts.footprint.RHO = 10;                % Density threshold
opts.footprint.PI = 0.75;               % Purity threshold
opts.footprint.LOWER_PCTILE = 1;        % Lower distance threshold
opts.footprint.UPPER_PCTILE = 25;       % Higher distance threshold

opts.pbldr.ntries = 10;                 % Number of attempts carried out by PBLDR
opts.pbldr.analytic = false;            % Calculate the analytical or numerical solution
opts.pbldr.cmaopts = bipopcmaes;        % Get the default params for BIPOP-CMA-ES
opts.pbldr.cmaopts.StopFitness = 0;     % Stop if the fitness is 0
opts.pbldr.cmaopts.MaxRestartFunEvals = 0;  % Allow multiple restarts
opts.pbldr.cmaopts.MaxFunEvals  = 1e4;      % Maximum number of evaluations
opts.pbldr.cmaopts.EvalParallel = 'no';     % This should be kept as 'no'
opts.pbldr.cmaopts.DispFinal = 'off';       % To make BIPOP-CMA-ES silent

% out = matilda(X, Y, Ybin, opts);
out = matilda(X(1:500,:), Y(1:500,:), Ybin(1:500,:), opts);

