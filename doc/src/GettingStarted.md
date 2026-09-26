# Getting Started with Instance Space Analysis

Build and explore your first instance space

This example builds an instance space for the reference dataset shipped with the toolbox: 212 classification problems (the *instances*), described by 10 features, on which 10 classification algorithms were scored by misclassification error (Muñoz et al., 2018). A further 23 instances are held out for testing. The same steps apply to your own data once it is in the [metadata format](MetadataFormat.html).

`example.m` in the repository root runs this example from start to finish.

## Prepare a Folder

The toolbox reads its inputs from, and writes its outputs to, one folder. Copy the reference metadata into a new folder:

```matlab
rootdir = 'test/data/example/';
if ~isfolder(rootdir), mkdir(rootdir); end
copyfile('test/data/metadata.csv', rootdir);
copyfile('test/data/metadata_test.csv', rootdir);
```

## Describe the Performance Measure

The one decision you must always make is how to judge performance. Here performance is an error rate, so lower is better, and an algorithm counts as *good* on an instance when its error is below 20%:

```matlab
opts.perf.MaxPerf = false;   % a cost: lower is better
opts.perf.AbsPerf = true;    % compare with a fixed threshold
opts.perf.epsilon = 0.20;    % the threshold
```

With `AbsPerf = false`, an algorithm is good when it is within `epsilon` (a fraction) of the best algorithm on that instance. Every other option has a default; see the [Options Reference](OptionsReference.html).

## Build the Instance Space

```matlab
obj = InstanceSpace(rootdir, opts);
obj = obj.build();
```

`build` runs the whole pipeline: `PRELIM` labels and normalises the data, `SIFTED` selects features, `PILOT` projects the instances to 2D, `CLOISTER` estimates the boundary of the space, `PYTHIA` trains a classifier per algorithm, and `TRACE` finds the footprints. It then saves `model.mat` and writes CSV and PNG files to `rootdir`. It takes a few minutes, most of it in SIFTED and PYTHIA.

## Look at the Results

```matlab
figure, obj.plot('portfolio')        % best algorithm on each instance
figure, obj.plot('footprint', 6)     % footprint of the 6th algorithm
disp(obj.model.trace.summary)        % footprint area, density and purity
disp(obj.model.pythia.summary)       % classifier accuracy and the selector
```

The PNG files in `rootdir` show the same views for every feature and algorithm.

## Evaluate New Instances

```matlab
obj = obj.explore(rootdir);          % reads rootdir/metadata_test.csv
test = obj.getResults(1);
scatter(test.pilot.Z(:,1), test.pilot.Z(:,2), 30, test.pythia.selection0, 'filled')
```

`explore` applies the trained preprocessing, projection, classifiers and footprints to the new instances; nothing is refitted.

## Next Steps

- Change an option and re-run only the affected stages — see `InstanceSpace`.
- Run and inspect each stage in turn — see the [Stage-by-Stage Walkthrough](InteractiveWalkthrough.html).
- Build a 3D instance space with `opts.pilot.dims = 3` — see `PILOT` and `PILOTviewpoint`.
- Use your own data — see [Metadata File Format](MetadataFormat.html).

## References

- Muñoz, M.A., Villanova, L., Baatar, D. & Smith-Miles, K. (2018). Instance spaces for machine learning classification. *Machine Learning*, 107(1), 109–147. <https://doi.org/10.1007/s10994-017-5629-5>
- Smith-Miles, K. & Muñoz, M.A. (2023). Instance Space Analysis for Algorithm Testing. *ACM Computing Surveys*, 55(12), Article 255. <https://doi.org/10.1145/3572895>
