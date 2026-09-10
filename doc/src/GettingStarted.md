# Getting Started

## Introduction

This page serves as a getting-started guide adapted from `example.m`, designed to help users familiarize themselves with the full Instance Space Analysis (ISA) pipeline. The tutorial covers the execution of `buildIS` and `exploreIS` using the bundled Munoz et al. 2018 reference dataset. This dataset consists of 212 instances, 10 features, and 10 classification algorithms, all scored by misclassification error.

## Running the Example

The following MATLAB code snippet from `example.m` sets up the environment and runs the ISA pipeline. Adjust the settings as necessary to suit your specific needs.

```matlab
srcdir  = './test/data/';
rootdir = './test/data/example/';
if ~isfolder(rootdir), mkdir(rootdir); end
copyfile([srcdir 'metadata.csv'], [rootdir 'metadata.csv']);
copyfile([srcdir 'metadata_test.csv'], [rootdir 'metadata_test.csv']);

% ---- A few of the most commonly adjusted settings --------------------------
classifier = 'knn';    % opts.pythia.classifier: 'knn' (default), 'svm', 'tree', 'nb', 'linear', 'ensemble'
tuning     = 'sobol';  % opts.pythia.tuning: 'sobol' (default), 'bayes', or 'none' (needs opts.pythia.params)
dims       = 2;        % opts.pilot.dims: 2 (default) or 3
siftedFlag = true;     % opts.sifted.flag: automated feature selection on/off
% -------------------------------------------------------------------------

opts = struct();
opts.pythia.classifier = classifier;
opts.pythia.tuning     = tuning;
opts.pilot.dims        = dims;
opts.sifted.flag       = siftedFlag;

% This dataset's algorithm performance is a misclassification error rate:
% lower is better (MaxPerf=false), and "good" means an absolute error
% below 20% (AbsPerf=true, epsilon=0.20). Everything else is left at the
% toolkit's defaults.
opts.perf.MaxPerf = false;
opts.perf.AbsPerf = true;
opts.perf.epsilon = 0.20;

fid = fopen([rootdir 'options.json'], 'w+');
fprintf(fid, '%s', jsonencode(opts));
fclose(fid);

model = buildIS(rootdir);
out = exploreIS(rootdir);
```

### Key Settings

- **`classifier`**: The classifier type, defaulting to 'knn'. Other options include 'svm', 'tree', 'nb', 'linear', and 'ensemble'.
- **`tuning`**: The tuning strategy, with options 'sobol' (default), 'bayes', or 'none' (requires additional parameters).
- **`dims`**: The projection dimensionality, defaulting to 2. The other option is 3.
- **`siftedFlag`**: A flag for automated feature selection, set to `true` for on.

### Performance Metric Settings

For this dataset, the performance metric is misclassification error. The goal is to minimize this error, where "good" performance is defined as an absolute error below 20%.

## Analyzing Custom Data

To analyze your own data, replace `SRCDIR` with the folder containing `metadata.csv`. Be sure to revisit the performance metric settings, as they are tuned for the reference dataset's error-rate semantics.

## Outputs

After running the script, the following workspace variables are available:

- **`model`**: The result of `buildIS`.
- **`out`**: The result of `exploreIS`.

Additionally, the script writes the following outputs to `rootdir`:

- CSVs
- PNGs
- `model.mat`

## See Also

For the exhaustive option-coverage regression suite, see `test_integration.m`. For the full list of configurable options, see the [Options Reference page](OptionsReference.html).

## Conclusion

This guide provides a step-by-step introduction to running the ISA pipeline using the bundled reference dataset. Adjust the settings as needed to fit your specific analysis requirements. For more advanced usage, refer to the regression suite and the full list of configurable options.
