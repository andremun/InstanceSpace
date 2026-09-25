# Stage-by-Stage Walkthrough

Run each stage of the pipeline and inspect its output

This walkthrough builds the instance space one stage at a time so you can see what each stage produces. It follows `liveDemoIS.m` in the repository root, which you can open in the Live Editor and run section by section. It uses the same reference data and performance settings as [Getting Started](GettingStarted.html).

```matlab
rootdir = 'test/data/liveDemo/';
if ~isfolder(rootdir), mkdir(rootdir); end
copyfile('test/data/metadata.csv', rootdir);
copyfile('test/data/metadata_test.csv', rootdir);

opts.perf = struct('MaxPerf', false, 'AbsPerf', true, 'epsilon', 0.20);
obj = InstanceSpace(rootdir, opts);
obj.completedStages          % empty: nothing has run yet
```

## PRELIM: Label and Normalise

`PRELIM` decides which algorithms are good on each instance, bounds outlying feature values, and normalises features and performance.

```matlab
obj = obj.build('stages', {'prelim'});
fprintf('%d instances, %d features, %d algorithms\n', ...
    size(obj.model.data.X,1), size(obj.model.data.X,2), size(obj.model.data.Y,2));
fprintf('Fraction of beta-easy instances: %.2f\n', mean(obj.model.data.beta));
```

An instance is *beta-easy* when more than `opts.perf.betaThreshold` of the algorithms are good on it.

## SIFTED: Select Features

`SIFTED` keeps the features that best explain performance: a correlation filter, then a genetic algorithm that picks one feature from each cluster of correlated features.

```matlab
obj = obj.build('stages', {'sifted'});
disp(obj.model.featsel.labels)       % before
disp(obj.model.data.featlabels)      % after
```

Set `opts.sifted.flag = false` to keep every feature.

## PILOT: Project to 2D

`PILOT` finds the linear projection `Z = X*A'` from which features and performance can best be reconstructed. `Z` is the instance space.

```matlab
obj = obj.build('stages', {'pilot'});
figure
scatter(obj.model.pilot.Z(:,1), obj.model.pilot.Z(:,2), 20, obj.model.data.numGoodAlgos, 'filled')
xlabel('z_1'), ylabel('z_2'), colorbar
title('Number of good algorithms per instance')
```

`opts.pilot.dims = 3` gives a 3D space; `opts.pilot.method = 'pls'` uses Partial Least Squares.

## CLOISTER: Estimate the Boundary

`CLOISTER` estimates where instances could exist, given the ranges and correlations of the features. Empty areas inside the boundary are candidates for new test instances.

```matlab
obj = obj.build('stages', {'cloister'});
figure, obj.plot('boundary')
```

## PYTHIA: Predict Good Algorithms

`PYTHIA` trains one classifier per algorithm that predicts, from an instance's position, whether the algorithm is good there, and combines them into an algorithm selector.

```matlab
obj = obj.build('stages', {'pythia'});
disp(obj.model.pythia.summary)
figure
scatter(obj.model.pilot.Z(:,1), obj.model.pilot.Z(:,2), 20, obj.model.pythia.selection0, 'filled')
colorbar, title('Selected algorithm')
```

## TRACE: Find the Footprints

`TRACE` finds the region where each algorithm is good (its footprint) and where it is best, and measures their size, density and purity.

```matlab
obj = obj.build('stages', {'trace'});
disp(obj.model.trace.summary)
figure, obj.plot('footprint', 6)
```

## Save and Explore

Once every stage has run, `build` has already saved `model.mat` and written the output files. `explore` evaluates new instances with the trained model:

```matlab
obj = obj.explore(rootdir);
test = obj.getResults(1);
disp(test.pythia.summary)
disp(test.trace.summary)
```

To change a setting later, edit `obj.opts` and re-run from the stage it affects, for example `obj = obj.build('stages', {'pythia','trace'})` after changing `opts.pythia.classifier`. Re-running a stage discards the results of later stages.

## See Also

`InstanceSpace` | [Getting Started](GettingStarted.html) | [Options Reference](OptionsReference.html)
