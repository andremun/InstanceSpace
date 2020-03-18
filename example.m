rootdir = 'E:/InstanceSpace_Classification/MATILDA_trial/';

opts.perf.MaxPerf = false;              % True if Y is a performance measure to maximize, False if it is a cost measure to minimise.
opts.perf.AbsPerf = true;               % True if an absolute performance measure, False if a relative performance measure
opts.perf.epsilon = 0.20;               % Threshold of good performance

opts.general.betaThreshold = 0.55;      % Beta-easy threshold

opts.auto.preproc = true;               % Automatic preprocessing on. Set to false if you don't want any preprocessing
opts.bound.flag = true;                 % Bound the outliers. True if you want to bound the outliers, false if you don't
opts.norm.flag = true;                  % Normalize/Standarize the data. True if you want to apply Box-Cox and Z transformations to stabilize the variance and scale N(0,1)

opts.auto.featsel = true;               % Automatic feature selectio on. Set to false if you don't want any feature selection.
opts.corr.flag = true;                  % Run feature selection by correlation between performance and features (Step 2)
opts.corr.threshold = 10;                % Top N features (by correlation) per algorithm that are selected
opts.clust.flag = true;                 % Run feature selection by clustering (Step 3)
opts.clust.KDEFAULT = 10;               % Default maximum number of clusters
opts.clust.SILTHRESHOLD = 0.9;          % Minimum accepted value for the average silhoute value
opts.clust.NTREES = 50;                 % Number of trees for the Random Forest (to determine highest separability in the 2-d projection)
opts.clust.MaxIter = 1000;
opts.clust.Replicates = 100;

opts.pilot.analytic = false;            % Calculate the analytical or numerical solution
opts.pilot.ntries = 5;                 % Number of attempts carried out by PBLDR

opts.cloister.pval = 0.05;
opts.cloister.cthres = 0.7;

opts.pythia.cvfolds = 5;
opts.pythia.useweights = false;

opts.trace.usesim = true;          % Use the actual or simulated data to calculate the footprints
opts.trace.RHO = 5;                 % Density threshold
opts.trace.PI = 0.70;               % Purity threshold

opts.selvars.smallscaleflag = false;     % True if you want to do a small scale experiment with a percentage of the available instances
opts.selvars.smallscale = 0.30;          % Percentage of instances to be kept for a small scale experiment
% You can also provide a file with the indexes of the instances to be used.
% This should be a csvfile with a single column of integer numbers that
% should be lower than the number of instances
opts.selvars.fileidxflag = false;
opts.selvars.fileidx = '';

opts.outputs.csv = true;               %
opts.outputs.web = true;              % NOTE: MAKE THIS FALSE IF YOU ARE USING THIS CODE LOCALY - This flag is only useful if the system is being used 'online' through matilda.unimelb.edu.au
opts.outputs.png = true;               %

% Saving all the information as a JSON file
fid = fopen([rootdir 'options.json'],'w+');
fprintf(fid,'%s',jsonencode(opts));
fclose(fid);

try
    model = trainIS(rootdir);
    testDir = strcat(rootdir, 'test/');
    testIS(testDir);
catch ME
    disp('EOF:ERROR');
    rethrow(ME)
end