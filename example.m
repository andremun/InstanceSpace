opts.perf.MaxMin = true;                % True if Y is a performance measure, False if it is a cost measure.
opts.perf.AbsPerf = true;               % True if an absolute performance measure, False if a relative performance measure
opts.perf.epsilon = 0.80;               % Threshold of good performance

opts.general.betaThreshold = 0.5;       % Beta-easy threshold

opts.bound.flag = true;

opts.norm.flag = true;

opts.diversity.flag = true;             % Run diversity calculation
opts.diversity.threshold = 0.3;         % Minimum percentage allowed of repeated values [?]

opts.corr.flag = true;
opts.corr.threshold = 3;                % Top N features (by correlation) per algorithm that are selected

opts.clust.flag = true;
opts.clust.KDEFAULT = 10;               % Default maximum number of clusters
opts.clust.SILTHRESHOLD = 0.58;          % Minimum accepted value for the average silhoute value
opts.clust.NTREES = 50;                 % Number of trees for the Random Forest (to determine highest separability in the 2-d projection)
opts.clust.MaxIter = 1000;
opts.clust.Replicates = 100;
opts.clust.UseParallel = true;

opts.oracle.cvgrid = 10;
opts.oracle.maxcvgrid = 5;
opts.oracle.mincvgrid = -5;
opts.oracle.cvfolds = 5;
opts.pbldr.analytic = false;            % Calculate the analytical or numerical solution
opts.pbldr.ntries = 30;                 % Number of attempts carried out by PBLDR

opts.footprint.RHO = 10;                % Density threshold
opts.footprint.PI = 0.75;               % Purity threshold
opts.footprint.LOWER_PCTILE = 5;        % Lower distance threshold
opts.footprint.UPPER_PCTILE = 25;       % Higher distance threshold

rootdir = 'E:\InstanceSpace_AnomalyDetection\ROC_Subset_Ftrs_Perfs\';
fid = fopen([rootdir 'options.json'],'w+');
fprintf(fid,'%s',jsonencode(opts));
fclose(fid);

try
    out = matilda([rootdir 'Features.csv'], ... % X
                  [rootdir 'perf_roc_subset.csv'],... % Y
                  [rootdir 'options.json'],... % options file
                  [rootdir 'ROC_Subset_Indices_3_Methods_Only.csv']); % can leave empty
catch ME
    disp('EOF:ERROR')
    rethrow(ME)
end