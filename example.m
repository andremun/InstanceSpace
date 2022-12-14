rootdir = './trial/';

opts.parallel.flag = false;
opts.parallel.ncores = 2;

opts.perf.MaxPerf = false;          % True if Y is a performance measure to maximize, False if it is a cost measure to minimise.
opts.perf.AbsPerf = true;           % True if an absolute performance measure, False if a relative performance measure
opts.perf.epsilon = 0.20;           % Threshold of good performance
opts.perf.betaThreshold = 0.55;     % Beta-easy threshold
opts.auto.preproc = true;           % Automatic preprocessing on. Set to false if you don't want any preprocessing
opts.bound.flag = true;             % Bound the outliers. True if you want to bound the outliers, false if you don't
opts.norm.flag = true;              % Normalize/Standarize the data. True if you want to apply Box-Cox and Z transformations to stabilize the variance and scale N(0,1)

opts.selvars.smallscaleflag = false; % True if you want to do a small scale experiment with a percentage of the available instances
opts.selvars.smallscale = 0.50;     % Percentage of instances to be kept for a small scale experiment
% You can also provide a file with the indexes of the instances to be used.
% This should be a csvfile with a single column of integer numbers that
% should be lower than the number of instances
opts.selvars.fileidxflag = false;
opts.selvars.fileidx = '';
opts.selvars.densityflag = false;
opts.selvars.mindistance = 0.1;
opts.selvars.type = 'Ftr&Good';

opts.sifted.flag = true;            % Automatic feature selectio on. Set to false if you don't want any feature selection.
opts.sifted.rho = 0.1;              % Minimum correlation value acceptable between performance and a feature. Between 0 and 1
opts.sifted.K = 10;                 % Number of final features. Ideally less than 10.
opts.sifted.NTREES = 50;            % Number of trees for the Random Forest (to determine highest separability in the 2-d projection)
opts.sifted.MaxIter = 1000;
opts.sifted.Replicates = 100;

opts.pilot.analytic = false;        % Calculate the analytical or numerical solution
opts.pilot.ntries = 5;              % Number of attempts carried out by PBLDR

opts.cloister.pval = 0.05;
opts.cloister.cthres = 0.7;

opts.pythia.cvfolds = 5;
opts.pythia.ispolykrnl = false;
opts.pythia.useweights = false;
opts.pythia.uselibsvm = false;

opts.trace.usesim = true;           % Use the actual or simulated data to calculate the footprints
opts.trace.PI = 0.55;               % Purity threshold

opts.outputs.csv = true;            %
opts.outputs.web = false;            % NOTE: MAKE THIS FALSE IF YOU ARE USING THIS CODE LOCALY - This flag is only useful if the system is being used 'online' through matilda.unimelb.edu.au
opts.outputs.png = true;            %

% Saving all the information as a JSON file
fid = fopen([rootdir 'options.json'],'w+');
fprintf(fid,'%s',jsonencode(opts));
fclose(fid);

try
    model = buildIS(rootdir);
catch ME
    disp('EOF:ERROR');
    rethrow(ME)
end