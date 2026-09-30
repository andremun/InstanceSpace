function isa_review_extra
addpath('/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace'); startup;
checks={@cacheCheck,@constantNorm,@cameraCheck,@seedCheck,@densityString,@tracePurity,@profilePilot};
for i=1:numel(checks)
 fprintf('\nREVIEW_CASE %s\n',func2str(checks{i}));
 try, checks{i}(); catch ME, fprintf('REVIEW_ERROR %s: %s\n',ME.identifier,ME.message); end
end
end
function cacheCheck
 rng(5); X=rand(40,6); Y=rand(40,2); Ybin=Y<.5;
 o=ISAdefaults(struct()); o.sifted.rho=0; o.sifted.K=3; o.sifted.Replicates=1;
 labels={'a','b','c','d','e','f'};
 clear SIFTED;
 profile clear; profile on; SIFTED(X,Y,Ybin,labels,o.sifted); profile off; a=profile('info');
 profile clear; profile on; SIFTED(X,Y,true(size(Ybin)),labels,o.sifted); profile off; b=profile('info');
 printCounts(a,'first'); printCounts(b,'second');
end
function printCounts(p,label)
 f=p.FunctionTable; ix=find(strcmp({f.FunctionName},'PILOT'));
 if isempty(ix), n=0; else,n=sum([f(ix).NumCalls]);end
 fprintf('REVIEW_RESULT %s SIFTED PILOT_calls=%d\n',label,n);
end
function constantNorm
 o=struct('MaxPerf',false,'AbsPerf',true,'epsilon',1,'betaThreshold',.55,'auto',true,'bound',true,'norm',true);
 X=[ones(20,1),rand(20,2)]; Y=rand(20,2);
 [a,b,p]=PRELIM(X,Y,o); [c,d]=PRELIM(X,Y,o,p);
 fprintf('REVIEW_RESULT constant feature train finite=%d eval finite=%d train_eval_X_diff=%g Y_diff=%g sigma=%s\n',all(isfinite(a),'all'),all(isfinite(c),'all'),max(abs(a-c),[],'all'),max(abs(b-d),[],'all'),mat2str(p.sigmaX));
end
function cameraCheck
 f=figure('Visible','off'); guard=onCleanup(@()close(f)); ax=axes(f); plot3(ax,[0 1],[0 1],[0 1]); axis(ax,'equal');
 targetDir=[1 0 0]; [az,el]=cart2sph(targetDir(1),targetDir(2),targetDir(3)); view(ax,rad2deg([az el]));
 actual=ax.CameraPosition-ax.CameraTarget; actual=actual/norm(actual);
 fprintf('REVIEW_RESULT desired direction=%s actual=%s\n',mat2str(targetDir),mat2str(actual,4));
end
function seedCheck
 o=ISAdefaults(struct()); o.pythia.seed=50000; o.pythia.nTuningIter=1; o.pythia.verbose=false;
 Ybin=repmat([true;false],10,1); PYTHIA(rand(20,2),rand(20,1),Ybin,rand(20,1),{'a'},o.pythia);
end
function densityString
 o.selvars.densityflag=true; o.selvars.type="Ftr"; ISAvalidateOpts(o);
 fprintf('REVIEW_RESULT validated scalar string passes bydensity type gate=%d\n',ischar(o.selvars.type));
end
function tracePurity
 Z=[0 0;1 0;0 1;1 1;.5 .5]; Z=[Z;Z]; y=[true(5,1);false(5,1)];
 o=ISAdefaults(struct()); o.trace.PI=.9; o.trace.minAreaFrac=0;
 a=TRACE(Z,y,y,ones(10,1),true(10,1),{'a'},o.trace);
 fprintf('REVIEW_RESULT accepted footprint measure=%g purity=%g threshold=%g\n',a.good{1}.measure,a.good{1}.purity,o.trace.PI);
end
function profilePilot
 rng(8); X=rand(600,5);Y=rand(600,2);o=ISAdefaults(struct());o.pilot.analytic=true;o.pilot.verbose=false;
 profile clear; profile on; PILOT(X,Y,{'a','b','c','d','e'},o.pilot);profile off;p=profile('info');
 f=p.FunctionTable; ix=contains({f.FunctionName},'pdist');fprintf('REVIEW_RESULT analytic branch pdist calls=%d\n',sum([f(ix).NumCalls]));
end
