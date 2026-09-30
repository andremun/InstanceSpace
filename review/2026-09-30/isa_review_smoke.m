function isa_review_smoke
addpath('/home/unimelb.edu.au/mariom1/Documents/packages/InstanceSpace'); startup;
d=[tempname '/'];mkdir(d);copyfile('test/data/metadata.csv',[d 'metadata.csv']);copyfile('test/data/metadata_test.csv',[d 'metadata_test.csv']);
o=ISAdefaults(struct());o.general.verbose=false;o.pilot.analytic=true;o.sifted.flag=false;o.outputs.csv=false;o.outputs.png=false;o.perf.AbsPerf=true;o.perf.epsilon=.2;o.pythia.nTuningIter=2;
a=InstanceSpace(d,o).build();b=InstanceSpace.load(d);b=b.explore(d);
fprintf('REVIEW_SMOKE train_rows=%d test_rows=%d classifiers=%d finite_train_Z=%d finite_test_Z=%d\n',size(a.model.pilot.Z,1),size(b.testResults{1}.pilot.Z,1),numel(a.model.pythia.classifiers),all(isfinite(a.model.pilot.Z),'all'),all(isfinite(b.testResults{1}.pilot.Z),'all'));
end
