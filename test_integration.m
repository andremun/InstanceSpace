% test_integration.m - ISA Toolkit full-pipeline regression suite (spec
% Section 12.2: "Full pipeline on reference dataset; numerical regression").
%
% This is the exhaustive option-coverage suite, not a getting-started
% example -- for that, see example.m, which runs a single bare-bones
% configuration with a short guide to the handful of settings most users
% adjust first. This file exists to catch regressions across the option
% surface: new features and bug fixes should get a case here (see
% docs/tech/isa_refactor_plan_v1.7.pdf, "Each phase is complete when... any
% new feature has a test case").
%
% Runs buildIS() + exploreIS() across a set of option configurations
% against the reference dataset (Munoz et al. 2018 classification study,
% test/data/metadata.csv + metadata_test.csv), each in its own subdirectory
% under test/data/ so outputs (model.mat, workspace.mat, CSVs, PNGs) never
% overwrite another case's results.
%
% IMPORTANT: options.json is NOT the source of truth for these settings.
% It is a generated artifact: this script writes a fresh options.json into
% each case's subdirectory (test/data/<case_name>/options.json) on every
% run, from the opts struct built in MATLAB below. Editing options.json
% directly has no lasting effect -- the next run of this script silently
% overwrites it. To change what gets tested, edit defaultOpts() or the
% relevant testCases{...}.override function in this file.
%
% Each test case is independent and wrapped in its own try/catch so one
% failure doesn't abort the rest of the suite; a summary is printed at the
% end, followed by the EOF:SUCCESS / EOF:ERROR sentinel already used by
% buildIS.m/exploreIS.m (so any automation scraping for that string still
% works against this script).

rootdir  = './test/data/';
srcfiles = {'metadata.csv', 'metadata_test.csv'};

% Number of algorithms in the reference dataset, read from the CSV header
% rather than hardcoded, so the 'tuning_none_precalc' case's pre-supplied
% opts.pythia.params matrix is always sized correctly even if the dataset
% changes.
hdr = readtable([rootdir 'metadata.csv']);
nalgos = sum(strncmpi(hdr.Properties.VariableNames, 'algo_', 5));

% ---- Shared baseline options -------------------------------------------
% Every test case starts here and overrides only the fields relevant to
% what it's checking. Tuning/trial counts are kept small (well below the
% spec defaults) purely so the whole suite runs quickly as a smoke test;
% increase them for a real experiment.
baseOpts = defaultOpts();

% ---- Test cases ----------------------------------------------------------
% Each entry: name (also the output subdirectory), a one-line description
% of what it exercises, and an override function applied to baseOpts.
testCases = {
    struct( ...
        'name', 'baseline_knn_sobol_2d', ...
        'desc', ['Default pipeline: KNN classifier, Sobol tuning, TRACE3, ' ...
                 'SIFTED on, 2D PILOT. The happy path every other case is a variation of.'], ...
        'override', @(o) o)

    struct( ...
        'name', 'classifier_svm', ...
        'desc', 'PYTHIA classifier registry: SVM (fitcsvm) instead of the KNN default.', ...
        'override', @(o) setPythia(o, 'classifier', 'svm'))

    struct( ...
        'name', 'classifier_ensemble', ...
        'desc', ['PYTHIA classifier registry: ensemble (fitcensemble, Bag method) -- ' ...
                 'the two-hyperparameter non-KNN case (NumLearningCycles/MinLeafSize).'], ...
        'override', @(o) setPythia(setPythia(o, 'classifier', 'ensemble'), 'ensembleMethod', 'Bag'))

    struct( ...
        'name', 'tuning_bayes', ...
        'desc', ['opts.pythia.tuning=''bayes'': MATLAB bayesopt hyperparameter search ' ...
                 'instead of the Sobol default.'], ...
        'override', @(o) setPythia(o, 'tuning', 'bayes'))

    struct( ...
        'name', 'tuning_none_precalc', ...
        'desc', ['opts.pythia.tuning=''none'' with pre-supplied opts.pythia.params: ' ...
                 'skips hyperparameter search entirely and uses fixed values.'], ...
        'override', @(o) setPythiaPrecalc(o, nalgos))

    struct( ...
        'name', 'pythia_skip', ...
        'desc', ['opts.pythia.skip=true: PYTHIA bypasses classifier training; TRACE must ' ...
                 'fall back to the true-label-only footprint (Zu={yi=1}) with a warning.'], ...
        'override', @(o) setPythia(o, 'skip', true))

    struct( ...
        'name', 'trace_legacy', ...
        'desc', 'opts.trace.method=''legacy'': the pre-refactor DBSCAN + alpha-shape footprint algorithm.', ...
        'override', @(o) setTrace(o, 'method', 'legacy'))

    struct( ...
        'name', 'sifted_off', ...
        'desc', 'opts.sifted.flag=false: pipeline runs without automated feature selection.', ...
        'override', @(o) setSifted(o, 'flag', false))

    struct( ...
        'name', 'small_scale_subset', ...
        'desc', 'opts.selvars.smallscaleflag=true: validates the small-scale subsetting path.', ...
        'override', @(o) setSmallScale(o, 0.5))

    struct( ...
        'name', 'density_subset', ...
        'desc', 'opts.selvars.densityflag=true: validates the FILTER-based density subsetting path.', ...
        'override', @(o) setDensitySubset(o))

    struct( ...
        'name', 'cost_sensitive_weights', ...
        'desc', ['opts.pythia.useweights=true: cost-sensitive classification weighted by ' ...
                 '|ybest-y|, including the degenerate-weights fallback guard.'], ...
        'override', @(o) setPythia(o, 'useweights', true))

    struct( ...
        'name', 'pilot_3d_numeric', ...
        'desc', ['opts.pilot.dims=3, opts.pilot.analytic=false: the numerical (BFGS) 3D projection ' ...
                 'path. Also exercises PILOTviewpoint''s default global-viewpoint path (buildIS calls ' ...
                 'it automatically whenever dims=3) and scriptpng''s use of the resulting azimuth/' ...
                 'elevation for every 3D plot (opts.outputs.png=true by default).'], ...
        'override', @(o) setPilot(o, 'dims', 3))

    struct( ...
        'name', 'pilot_analytic_2d', ...
        'desc', 'opts.pilot.analytic=true, 2D: the closed-form eigenvector solution (top-2 eigenvectors).', ...
        'override', @(o) setPilot(o, 'analytic', true))

    struct( ...
        'name', 'pilot_analytic_3d', ...
        'desc', ['opts.pilot.analytic=true, opts.pilot.dims=3: the generalised closed-form ' ...
                 'solution (top-3 eigenvectors) -- exercises the 3D analytic branch specifically.'], ...
        'override', @(o) setPilot(setPilot(o, 'analytic', true), 'dims', 3))

    struct( ...
        'name', 'pilot_isa3d_legacy_alias', ...
        'desc', ['opts.pilot.ISA3D=true with dims removed: regression check that the legacy boolean ' ...
                 'alias still maps to the 3D path via ISAdefaults/PILOT''s backward-compat shim.'], ...
        'override', @(o) setPilotISA3DLegacy(o))

    struct( ...
        'name', 'pilot_pls_2d', ...
        'desc', 'opts.pilot.method=''pls'', 2D: the Partial Least Squares alternative projection method.', ...
        'override', @(o) setPilot(o, 'method', 'pls'))

    struct( ...
        'name', 'pilot_pls_3d', ...
        'desc', 'opts.pilot.method=''pls'', opts.pilot.dims=3: PLS at 3D.', ...
        'override', @(o) setPilot(setPilot(o, 'method', 'pls'), 'dims', 3))

    struct( ...
        'name', 'pilot_alpha_weight', ...
        'desc', ['opts.pilot.alpha=3.0: emphasises the performance-reconstruction term over the ' ...
                 'feature-reconstruction term in PILOT''s cost function.'], ...
        'override', @(o) setPilot(o, 'alpha', 3.0))

    struct( ...
        'name', 'pilot_viewpoint_groups', ...
        'desc', ['opts.pilot.dims=3, opts.pilot.viewGroups={1:4,5:10}: exercises PILOTviewpoint''s ' ...
                 'per-group path (one viewpoint per algorithm group instead of a single global one), ' ...
                 'and scriptpng''s per-algorithm group lookup (resolveViewAngle) that picks the right ' ...
                 'viewpoint for each algorithm''s plots. Uneven group sizes so the options.json ' ...
                 'round-trip decodes as a cell array rather than collapsing into a numeric matrix ' ...
                 '(jsonencode/jsondecode quirk).'], ...
        'override', @(o) setPilot(setPilot(o, 'dims', 3), 'viewGroups', {1:4, 5:10}))

    struct( ...
        'name', 'feats_algos_subset', ...
        'desc', ['opts.selvars.feats/algos: restrict to a hand-picked subset of features ' ...
                 'and algorithms (README-documented, must be preserved).'], ...
        'override', @(o) setFeatsAlgosSubset(o))

    struct( ...
        'name', 'sifted_parallel_ga', ...
        'desc', ['opts.general.parallel=true: exercises SIFTED''s GA-level UseParallel path ' ...
                 'and clearCache''s cross-worker persistent-cache reset via parfevalOnAll ' ...
                 '(spec Phase 6). Every other case leaves opts.general.parallel=false, so ' ...
                 'this path was previously untested. Kept last in the list since it opens a ' ...
                 'parallel pool that stays open for the remainder of the script.'], ...
        'override', @(o) setGeneralParallel(o))
};

% ---- Run all cases --------------------------------------------------------
nCases  = numel(testCases);
results = struct('name', cell(1,nCases), 'passed', cell(1,nCases), 'message', cell(1,nCases));

for i = 1:nCases
    tc = testCases{i};
    fprintf('\n[EXAMPLE] === Case %d/%d: %s ===\n', i, nCases, tc.name);
    fprintf('[EXAMPLE] %s\n', tc.desc);

    caseDir = [rootdir tc.name '/'];
    if ~isfolder(caseDir), mkdir(caseDir); end
    for f = 1:numel(srcfiles)
        copyfile([rootdir srcfiles{f}], [caseDir srcfiles{f}]);
    end

    opts = tc.override(baseOpts);
    fid = fopen([caseDir 'options.json'], 'w+');
    fprintf(fid, '%s', jsonencode(opts));
    fclose(fid);

    results(i).name = tc.name;
    try
        model = buildIS(caseDir); %#ok<NASGU>
        out = exploreIS(caseDir); %#ok<NASGU>
        results(i).passed  = true;
        results(i).message = 'OK';
        fprintf('[EXAMPLE] Case ''%s'' PASSED.\n', tc.name);
    catch ME
        results(i).passed  = false;
        results(i).message = ME.message;
        fprintf('[EXAMPLE] Case ''%s'' FAILED: %s\n', tc.name, ME.message);
    end
end

% ---- Summary ----------------------------------------------------------
fprintf('\n[EXAMPLE] ================= Summary =================\n');
nPassed = 0;
for i = 1:nCases
    if results(i).passed
        status = 'PASS';
        nPassed = nPassed + 1;
    else
        status = 'FAIL';
    end
    fprintf('[EXAMPLE] [%s] %s\n', status, results(i).name);
end
fprintf('[EXAMPLE] %d/%d cases passed.\n', nPassed, nCases);

if nPassed == nCases
    fprintf('EOF:SUCCESS\n');
else
    fprintf('EOF:ERROR\n');
    error('ISA:test_integration:caseFailures', '%d of %d test cases failed.', nCases-nPassed, nCases);
end

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function opts = defaultOpts()
% Baseline options shared by every test case below. Kept close to the
% published defaults (see README.md / docs/tech/isa_refactor_plan_v1.7.pdf
% Appendix A) except where noted.
opts.general.parallel = false;
opts.general.ncores = 18;

opts.perf.MaxPerf = false;          % True if Y is a performance measure to maximize, False if it is a cost measure to minimise.
opts.perf.AbsPerf = true;           % True if an absolute performance measure, False if a relative performance measure
opts.perf.epsilon = 0.20;           % Threshold of good performance
opts.perf.betaThreshold = 0.55;     % Beta-easy threshold
opts.auto.preproc = true;           % Automatic preprocessing on. Set to false if you don't want any preprocessing
opts.bound.flag = true;             % Bound the outliers. True if you want to bound the outliers, false if you don't
opts.norm.flag = true;              % Normalize/Standarize the data. True if you want to apply Box-Cox and Z transformations to stabilize the variance and scale N(0,1)

opts.selvars.smallscaleflag = false;
opts.selvars.smallscale = 0.50;
opts.selvars.fileidxflag = false;
opts.selvars.fileidx = '';
opts.selvars.densityflag = false;
opts.selvars.mindistance = 0.1;
opts.selvars.type = 'Ftr&Good';

opts.sifted.flag = true;
opts.sifted.rho = 0.1;              % Minimum correlation value acceptable between performance and a feature. Between 0 and 1
% K=5 (not the default 10) so SIFTED's k-means+GA clustering path is
% actually exercised: the reference dataset has only 10 features, and
% K=nfeats triggers SIFTED's "fewer features than clusters" skip branch
% instead (spec §12.1 testing guidance: "SIFTED may be tested with
% opts.sifted.K = 5 or 6").
opts.sifted.K = 5;
opts.sifted.MaxIter = 1000;
opts.sifted.Replicates = 100;

opts.pilot.analytic = false;        % Numerical (BFGS) by default; specific cases below flip this on.
opts.pilot.ntries = 3;              % Small trial count for a quick smoke test; increase for a real run.
opts.pilot.dims = 2;                % 2D by default; specific cases below flip this to 3.
opts.pilot.method = 'standard';     % 'standard' (BFGS/analytic) or 'pls'; specific cases below flip this.
opts.pilot.alpha = 1.0;             % Performance-reconstruction cost weight; specific cases below vary this.

opts.cloister.pval = 0.05;
opts.cloister.corrThreshold = 0.7;

opts.pythia.flag = true;
opts.pythia.classifier = 'knn';      % Registry name: 'knn','svm','tree','nb','linear','ensemble'
opts.pythia.tuning = 'sobol';
opts.pythia.nTuningIter = 5;         % Small evaluation budget for a quick smoke test; increase for a real run.
opts.pythia.kFold = 5;
opts.pythia.ispolykrnl = false;
opts.pythia.useweights = false;
opts.pythia.skip = false;

opts.trace.method = 'trace3';
opts.trace.PI = 0.55;               % Purity threshold

opts.outputs.csv = true;
opts.outputs.web = false;           % NOTE: MAKE THIS FALSE IF YOU ARE USING THIS CODE LOCALLY - This flag is only useful if the system is being used 'online' through matilda.unimelb.edu.au
opts.outputs.png = true;
end

function opts = setPythia(opts, field, value)
opts.pythia.(field) = value;
end

function opts = setSifted(opts, field, value)
opts.sifted.(field) = value;
end

function opts = setPilot(opts, field, value)
opts.pilot.(field) = value;
end

function opts = setGeneralParallel(opts)
opts.general.parallel = true;
opts.general.ncores = 2; % small pool, just enough to exercise the parallel path
end

function opts = setPilotISA3DLegacy(opts)
% Removes the dims field (set unconditionally in defaultOpts) so
% ISAdefaults' legacy-alias migration -- which only fires when dims is
% absent -- actually gets exercised by this case.
opts.pilot = rmfield(opts.pilot, 'dims');
opts.pilot.ISA3D = true;
end

function opts = setTrace(opts, field, value)
opts.trace.(field) = value;
end

function opts = setPythiaPrecalc(opts, nalgos)
opts.pythia.tuning = 'none';
% Fixed KNN hyperparameters for every algorithm: NumNeighbors=5, Distance
% index 1 (euclidean, per ISAgetClassifierFcn's distance ordering).
% nalgos is the raw algorithm count from metadata.csv's header (see top of
% this script); this assumes buildIS.m's "algorithms with no good instances"
% filter doesn't drop any of them under this case's perf/epsilon settings.
% If it ever does, PYTHIA's own validation raises a clear
% ISA:PYTHIA:noParamsForNoneTuning error rather than failing silently.
opts.pythia.params = [repmat(5, nalgos, 1), repmat(1, nalgos, 1)];
end

function opts = setSmallScale(opts, fraction)
opts.selvars.smallscaleflag = true;
opts.selvars.smallscale = fraction;
end

function opts = setDensitySubset(opts)
opts.selvars.densityflag = true;
opts.selvars.mindistance = 0.1;
opts.selvars.type = 'Ftr&Good';
end

function opts = setFeatsAlgosSubset(opts)
% Names must match the raw CSV header (feature_/algo_ prefixes intact):
% buildIS.m applies this filter before stripping the prefixes.
opts.selvars.feats = {'feature_Max_Normalized_Entropy_attributes', ...
                       'feature_ErrorRate_Decision_Node', ...
                       'feature_Max_Feature_Efficiency_F3', ...
                       'feature_Training_Error_Linear_Classifier_L2', ...
                       'feature_Nonlinearity_Nearest_Neighbor_Classifier_N4'};
opts.selvars.algos = {'algo_NB', 'algo_CART', 'algo_KNN'};
end
