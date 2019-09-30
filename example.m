rootdir = './';

opts.perf.MaxPerf = false;             % True if Y is a performance measure to maximize, False if it is a cost measure to minimise.
opts.perf.AbsPerf = false;              % True if an absolute performance measure, False if a relative performance measure
opts.perf.epsilon = 0.05;               % Threshold of good performance

opts.general.betaThreshold = 0.55;      % Beta-easy threshold

opts.auto.preproc = true;               % Automatic preprocessing on. Set to false if you don't want any preprocessing
opts.bound.flag = true;                 % Bound the outliers. True if you want to bound the outliers, false if you don't
opts.norm.flag = true;                  % Normalize/Standarize the data. True if you want to apply Box-Cox and Z transformations to stabilize the variance and scale N(0,1)

opts.auto.featsel = true;               % Automatic feature selectio on. Set to false if you don't want any feature selection.
opts.diversity.flag = false;            % Run feature selection by diversity calculation (Step 1)
opts.diversity.threshold = 0.3;         % Minimum percentage allowed of repeated values before a feature is considered NOT diverse enough
opts.corr.flag = true;                  % Run feature selection by correlation between performance and features (Step 2)
opts.corr.threshold = 6;                % Top N features (by correlation) per algorithm that are selected
opts.clust.flag = true;                 % Run feature selection by clustering (Step 3)
opts.clust.KDEFAULT = 10;               % Default maximum number of clusters
opts.clust.SILTHRESHOLD = 0.9;          % Minimum accepted value for the average silhoute value
opts.clust.NTREES = 50;                 % Number of trees for the Random Forest (to determine highest separability in the 2-d projection)
opts.clust.MaxIter = 1000;
opts.clust.Replicates = 100;

opts.pbldr.analytic = false;            % Calculate the analytical or numerical solution
opts.pbldr.ntries = 10;                 % Number of attempts carried out by PBLDR

opts.sbound.pval = 0.05;
opts.sbound.cthres = 0.7;

opts.oracle.cvgrid = 10;
opts.oracle.maxcvgrid = 5;
opts.oracle.mincvgrid = -5;
opts.oracle.cvfolds = 5;

opts.footprint.usesim = false;          % Use the actual or simulated data to calculate the footprints
opts.footprint.RHO = 10;                % Density threshold
opts.footprint.PI = 0.75;               % Purity threshold
opts.footprint.PCTILE = 0.30;

opts.selvars.smallscaleflag = true;      % True if you want to do a small scale experiment with a percentage of the available instances
opts.selvars.smallscale = 0.30;          % Percentage of instances to be kept for a small scale experiment
% You can also provide a file with the indexes of the instances to be used.
% This should be a csvfile with a single column of integer numbers that
% should be lower than the number of instances
opts.selvars.fileidxflag = false;
opts.selvars.fileidx = '';

opts.outputs.csv = true;               %
opts.outputs.web = false;              % NOTE: MAKE THIS FALSE IF YOU ARE USING THIS CODE LOCALY - This flag is only useful if the system is being used 'online' through matilda.unimelb.edu.au
opts.outputs.png = true;               %

% Saving all the information as a JSON file
fid = fopen([rootdir 'options.json'],'w+');
fprintf(fid,'%s',jsonencode(opts));
fclose(fid);

try
    model = trainIS(rootdir);
    testIS(rootdir);
catch ME
    disp('EOF:ERROR');
    rethrow(ME)
end