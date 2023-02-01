rootdir = 'C:/Users/mariom1/OneDrive - The University of Melbourne/Documents/MATLAB/InstanceSpace_Regression/trial/';

optsfile = [rootdir 'options_bkup.json'];
opts = jsondecode(fileread(optsfile));
opts.auto.featsel = false;
opts.corr.flag = true;
opts.corr.threshold = 5;
opts.clust.flag = false;

opts.pilot = opts.pbldr;
opts.pilot.ntries = 30;
opts.pilot.alpha = [0.565534791861630;0.530308089016240;-0.158055472179465;0.441925817765667;-0.388131680786894;-0.115897084892352;0.318031215594259;-0.252996348043433;-0.413799852670937;-0.0986568416066856;0.288430949599239;-0.451900231482780;0.0920288306088186;0.210208754444068;0.501893783783401;-0.514585346888767;-0.511216943113022;0.423274119894000;-0.454752084097749;0.265861589776303;0.401368888352203;-0.310344678529976;-0.349567685663719;-0.236589491913245;-0.287934053324859;-0.395479208893403;-0.484583421210239;-0.492561449063099;-0.335035242672287;0.482304931836019;0.261015742474904;-0.146334495978856;-0.421128618719782;-0.0470575423262341;-0.808393669073297;0.644318607634159;-0.304705280866874;-0.441763092040982;-0.292545080896145;-0.325736819687666;-0.526130539808091;-0.368615647760456;-0.424904810984364;-0.411911985220851];
opts.cloister = opts.sbound;
opts.pythia = opts.oracle;
opts.pythia.useweights = false;
opts.trace = opts.footprint;
opts.trace.PI = 0.65;

opts.selvars.fileidxflag = true;
opts.selvars.fileidx = [rootdir 'indexes.csv'];
opts.selvars.algos = {'algo_grad_boost';
                      'algo_bayesian_ARD';
                      'algo_bagging';
                      'algo_random_forest';
                      'algo_nu_svr';
                      'algo_linear_svr';
                      'algo_adaboost';
                      'algo_decision_tree'};
                  
opts.selvars.feats = {'feature_n1';
                      'feature_c2';
                      'feature_c5';
                      'feature_l1_a';
                      'feature_m5';
                      'feature_s1';
                      'feature_t2'};

fid = fopen([rootdir 'options.json'],'w+');
fprintf(fid,'%s',jsonencode(opts));
fclose(fid);

try
    model = trainIS(rootdir);
    PI=0:0.05:0.95;
    nalgos = length(model.data.algolabels);
    summaryExp = zeros(nalgos,10,length(PI));
    summarySim = zeros(nalgos,10,length(PI));
    for i=1:length(PI)
        opts.trace.PI = PI(i);
        model.trace.exper{i} = TRACE(model.pilot.Z, model.data.Ybin, model.data.P, model.data.beta, model.data.algolabels, opts.trace);
        model.trace.simul{i} = TRACE(model.pilot.Z, model.pythia.Yhat, model.pythia.selection0, model.data.beta, model.data.algolabels, opts.trace);
        summaryExp(:,:,i) = cell2mat(model.trace.exper{i}.summary(2:end,2:end));
        summarySim(:,:,i) = cell2mat(model.trace.simul{i}.summary(2:end,2:end));
    end
    summary{1} = [PI' squeeze(mean(summaryExp(:,[2 4 5 7 9 10],:),1))'];
    summary{2} = [PI' squeeze(mean(summarySim(:,[2 4 5 7 9 10],:),1))'];
catch ME
    disp('EOF:ERROR');
    rethrow(ME)
end

