function isa_review_checks
addpath('/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace'); startup;
checks={@pruning,@partialSave,@plsProjection,@smallEval,@fallbackLeak,@recallMetric,@weights,@jsonGroups,@missingFeature,@siftedRerun,@nanPilot,@negativePerf,@rankFallback,@stringRoot};
for i=1:numel(checks)
 fprintf('\nREVIEW_CASE %s\n',func2str(checks{i}));
 try, checks{i}(); catch ME, fprintf('REVIEW_ERROR %s: %s\n',ME.identifier,ME.message); end
end
end
function [o,d]=fixture
 d=[tempname '/']; mkdir(d); rng(4); n=40;
 T=array2table([rand(n,4)+2, 5+rand(n,1), .1+rand(n,1)*.1, .01+rand(n,1)*.02], 'VariableNames',{'feature_a','feature_b','feature_c','feature_d','algo_drop','algo_b','algo_c'});
 T=addvars(T,cellstr(string((1:n)')),'Before',1,'NewVariableNames','instances');
 writetable(T,[d 'metadata.csv']); writetable(T,[d 'metadata_test.csv']);
 o=ISAdefaults(struct()); o.general.verbose=false; o.auto.preproc=false; o.sifted.flag=false; o.pilot.analytic=true;
 o.perf.AbsPerf=true; o.perf.epsilon=1; o.outputs.csv=false; o.outputs.png=false; o.pythia.skip=true;
end
function pruning
 [o,d]=fixture; a=InstanceSpace(d,o).build('stages',{'prelim'});
 fprintf('REVIEW_RESULT retained=%d P=%s beta_true=%d expected_beta_true=%d\n',numel(a.model.data.algolabels),mat2str(unique(a.model.data.P)'),sum(a.model.data.beta),sum(sum(a.model.data.Ybin,2)>o.perf.betaThreshold*size(a.model.data.Ybin,2)));
end
function partialSave
 [o,d]=fixture; a=InstanceSpace(d,o).build('stages',{'prelim'}); a.save(); b=InstanceSpace.load(d); disp(b.completedStages);
end
function plsProjection
 rng(1); X=rand(30,4)+10; Y=rand(30,2); o=ISAdefaults(struct()); o.pilot.method='pls'; o.pilot.verbose=false;
 p=PILOT(X,Y,{'a','b','c','d'},o.pilot);
 fprintf('REVIEW_RESULT PLS train/eval maximum difference %.9g\n',max(abs(p.Z-X*p.A'),[],'all'));
end
function smallEval
 rng(3); Z=rand(30,2); o=ISAdefaults(struct()); t=TRACE(Z,true(30,1),true(30,1),ones(30,1),true(30,1),{'a'},o.trace);
 TRACE(Z(1,:),true,true,1,true,{'a'},o.trace,t);
end
function fallbackLeak
 t.mu=[0 0]; t.sigma=[1 1]; t.classifiers={struct('constant',true,'value',false),struct('constant',true,'value',false)}; t.precision=[.8;.7];
 Z=rand(8,2); Y=rand(8,2); y1=repmat([true false],8,1); y2=~y1;
 a=PYTHIA(Z,Y,y1,min(Y,[],2),{'a','b'},struct(),t); b=PYTHIA(Z,Y,y2,min(Y,[],2),{'a','b'},struct(),t);
 fprintf('REVIEW_RESULT same features/model different truth selects %d vs %d\n',a.selection1(1),b.selection1(1));
end
function recallMetric
 Z=rand(10,2); Y=rand(10,2); o=ISAdefaults(struct()); a=PYTHIA(Z,Y,true(10,2),min(Y,[],2),{'a','b'},o.pythia);
 fprintf('REVIEW_RESULT all predictions and selections correct selector recall=%g\n',a.summary{end,9});
end
function weights
 Z=rand(10,2); Y=[(1:10)' (11:20)']; o=ISAdefaults(struct()); o.pythia.useweights=true;
 a=PYTHIA(Z,Y,true(10,2),min(Y,[],2),{'a','b'},o.pythia);
 fprintf('REVIEW_RESULT weights first row actual=%s regret=%s\n',mat2str(a.W(1,:)),mat2str(abs(Y(1,:)-min(Y(1,:)))));
end
function jsonGroups
 o.pilot.viewGroups={[1 2],[3 4]}; q=jsondecode(jsonencode(o)); disp(q.pilot.viewGroups); ISAvalidateOpts(q);
end
function missingFeature
 [o,d]=fixture; T=readtable([d 'metadata.csv']); T.feature_a(1:10)=NaN; writetable(T,[d 'metadata.csv']);
 a=InstanceSpace(d,o).build('stages',{'prelim'}); INIT(d,a.opts,a.model);
end
function siftedRerun
 [o,d]=fixture; o.sifted.flag=true; o.sifted.rho=1; o.sifted.K=10;
 a=InstanceSpace(d,o).build('stages',{'prelim','sifted'}); n1=size(a.model.data.X,2);
 a.opts.sifted.flag=false; a=a.build('stages',{'sifted'});
 fprintf('REVIEW_RESULT original 4 features after sifted=%d after disabling sifted=%d\n',n1,size(a.model.data.X,2));
end
function nanPilot
 o=ISAdefaults(struct()); o.pilot.analytic=true; X=rand(20,4); X(1,1)=NaN;
 PILOT(X,rand(20,2),{'a','b','c','d'},o.pilot);
end
function negativePerf
 o=struct('MaxPerf',false,'AbsPerf',false,'epsilon',.05,'betaThreshold',.55,'auto',false,'bound',false,'norm',false);
 [~,~,p]=PRELIM(rand(3,3),repmat([-10 -1],3,1),o); fprintf('REVIEW_RESULT negative min costs [-10 -1] Ybin=%s\n',mat2str(p.Ybin(1,:)));
 o.MaxPerf=true; [~,~,p]=PRELIM(rand(3,3),zeros(3,2),o); fprintf('REVIEW_RESULT maximize all zeros Ybin=%s\n',mat2str(p.Ybin(1,:)));
end
function rankFallback
 X=rand(20,2); X=[X X(:,1)]; PILOT(X,rand(20,2),{'a','b','c'},struct('analytic',true,'dims',2,'verbose',false));
end
function stringRoot
 [o,d]=fixture; a=InstanceSpace(string(d(1:end-1)),o); disp(a.rootdir);
end
