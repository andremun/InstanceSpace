% test_integration.m - ISA Toolkit full-pipeline regression suite: runs the
% full pipeline on the reference dataset and checks for numerical regressions.
%
% This is the exhaustive option-coverage suite, not a getting-started
% example -- for that, see example.m, which runs a single bare-bones
% configuration with a short guide to the handful of settings most users
% adjust first. This file exists to catch regressions across the option
% surface: new features and bug fixes should get a case here.
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

% -------------------------------------------------------------------------
% Instance Space Analysis (ISA) Toolkit
% Copyright (c) 2026 Mario Andres Munoz Acosta and contributors
% School of Computing and Information Systems
% The University of Melbourne, Australia
%
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
% License: https://polyformproject.org/licenses/noncommercial/1.0.0/
%
% You may use, modify, and redistribute this software for non-commercial
% research and educational purposes only. Commercial use requires prior
% written permission. See the LICENSE file for full terms.
%
% Reference:
%   Smith-Miles, K. & Munoz, M.A. (2023). Instance Space Analysis for
%   Algorithm Testing. ACM Computing Surveys, 55(12), Article 255.
%   https://doi.org/10.1145/3572895
% -------------------------------------------------------------------------

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
    fprintf('\n[TEST] === Case %d/%d: %s ===\n', i, nCases, tc.name);
    fprintf('[TEST] %s\n', tc.desc);

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
        fprintf('[TEST] Case ''%s'' PASSED.\n', tc.name);
    catch ME
        results(i).passed  = false;
        results(i).message = ME.message;
        fprintf('[TEST] Case ''%s'' FAILED: %s\n', tc.name, ME.message);
    end
end

% ---- InstanceSpace class API coverage (spec Phase 7) -----------------
% The cases above only exercise InstanceSpace indirectly, through the
% buildIS/exploreIS backward-compatibility wrappers (each of which always
% runs every stage in one shot). This drives the class directly: staged
% build() with an option change between stages, out-of-order stage
% requests, the missing-prerequisite error path, a save()/load()
% round-trip, and explore() without going through exploreIS.m.
fprintf('\n[TEST] === Class API: staged build + save/load + explore ===\n');
classCaseDir = [rootdir 'class_api/'];
if ~isfolder(classCaseDir), mkdir(classCaseDir); end
for f = 1:numel(srcfiles)
    copyfile([rootdir srcfiles{f}], [classCaseDir srcfiles{f}]);
end
results(end+1).name = 'class_api';
try
    obj = InstanceSpace(classCaseDir, baseOpts);

    obj = obj.build('stages', {'prelim', 'sifted', 'pilot'});
    assert(isequal(obj.completedStages, {'prelim', 'sifted', 'pilot'}), ...
        'completedStages mismatch after a partial build().');

    obj.opts.pilot.alpha = 2.0;
    obj = obj.build('stages', {'pilot'}); % re-run just PILOT with the new weight
    assert(isequal(obj.completedStages, {'prelim', 'sifted', 'pilot'}), ...
        're-running an already-completed stage should not duplicate it in completedStages.');

    % Requested out of canonical order: build() must still run them
    % prelim->...->trace internally regardless of the order listed here.
    obj = obj.build('stages', {'trace', 'cloister', 'pythia'});
    assert(all(ismember({'cloister', 'pythia', 'trace'}, obj.completedStages)), ...
        'cloister/pythia/trace did not complete.');

    % Re-running an EARLIER stage after later ones have already completed
    % must invalidate cloister/pythia/trace, not leave them looking valid
    % alongside a freshly re-run pilot.
    obj.opts.pilot.alpha = 3.0;
    obj = obj.build('stages', {'pilot'});
    assert(isequal(obj.completedStages, {'prelim', 'sifted', 'pilot'}), ...
        're-running pilot after cloister/pythia/trace completed should invalidate them in completedStages.');
    assert(~isfield(obj.model, 'cloist') && ~isfield(obj.model, 'pythia') && ~isfield(obj.model, 'trace'), ...
        're-running pilot should remove the now-stale cloister/pythia/trace model fields.');

    % explore() must refuse a model left partially invalidated like this,
    % rather than crash deep inside evaluateTestSet on a missing field.
    notBuiltEnforced = false;
    try
        obj.explore(classCaseDir);
    catch notBuiltErr
        notBuiltEnforced = strcmp(notBuiltErr.identifier, 'ISA:InstanceSpace:notBuilt');
    end
    assert(notBuiltEnforced, 'explore() on a partially-invalidated model should raise ISA:InstanceSpace:notBuilt.');

    % Complete the pipeline again before the save()/load()/explore() checks below.
    obj = obj.build('stages', {'cloister', 'pythia', 'trace'});

    % Re-running 'cloister' must NOT invalidate 'pythia'/'trace': per
    % StagePrereq neither depends on cloister's output (both depend on
    % 'pilot'/'pythia' instead), even though both appear later than
    % 'cloister' in canonical StageOrder.
    pythiaBefore = obj.model.pythia;
    traceBefore = obj.model.trace;
    obj = obj.build('stages', {'cloister'});
    assert(all(ismember({'cloister', 'pythia', 'trace'}, obj.completedStages)), ...
        're-running cloister should not invalidate pythia/trace in completedStages.');
    assert(isequal(obj.model.pythia, pythiaBefore) && isequal(obj.model.trace, traceBefore), ...
        're-running cloister should not recompute/discard the still-valid pythia/trace results.');

    % A genuinely missing prerequisite must error clearly, not crash deep
    % inside the requested stage.
    prereqEnforced = false;
    try
        InstanceSpace(classCaseDir, baseOpts).build('stages', {'pilot'});
    catch prereqErr
        prereqEnforced = strcmp(prereqErr.identifier, 'ISA:InstanceSpace:missingPrereq');
    end
    assert(prereqEnforced, 'build(''stages'',{''pilot''}) on a fresh object should raise ISA:InstanceSpace:missingPrereq.');

    % save()/load() round-trip.
    obj.save();
    loaded = InstanceSpace.load(classCaseDir);
    assert(isequal(loaded.model.pilot.A, obj.model.pilot.A), 'save()/load() round-trip changed model.pilot.A.');

    % explore() directly through the class, not via exploreIS.m.
    loaded = loaded.explore(classCaseDir);
    assert(numel(loaded.testResults) == 1 && strcmp(loaded.testDirs{1}, classCaseDir), ...
        'explore() did not record testResults/testDirs correctly.');
    assert(isfield(loaded.getResults(1), 'trace'), 'getResults(1) is missing the trace field.');

    results(end).passed  = true;
    results(end).message = 'OK';
    fprintf('[TEST] Case ''class_api'' PASSED.\n');
catch ME
    results(end).passed  = false;
    results(end).message = ME.message;
    fprintf('[TEST] Case ''class_api'' FAILED: %s\n', ME.message);
end

% ---- ISAmigrateModel legacy migration coverage (spec §6.4) -----------
fprintf('\n[TEST] === ISAmigrateModel: legacy migration table ===\n');
results(end+1).name = 'isamigratemodel';
try
    % Pure opts/data field migrations: synthetic legacy-shaped structs, no
    % trained model needed.
    m = ISAmigrateModel(struct('opts', struct( ...
        'pbldr',     struct('dims', 2), ...
        'sbound',    struct('pval', 0.05), ...
        'footprint', struct('method', 'trace3'))));
    assert(isfield(m.opts, 'pilot') && m.opts.pilot.dims == 2 && ~isfield(m.opts, 'pbldr'), ...
        'opts.pbldr -> opts.pilot rename failed.');
    assert(isfield(m.opts, 'cloister') && ~isfield(m.opts, 'sbound'), ...
        'opts.sbound -> opts.cloister rename failed.');
    assert(isfield(m.opts, 'trace') && ~isfield(m.opts, 'footprint'), ...
        'opts.footprint -> opts.trace rename failed.');

    m = ISAmigrateModel(struct('opts', struct( ...
        'corr',  struct('flag', true, 'threshold', 0.35), ...
        'clust', struct('flag', false))));
    assert(~isfield(m.opts, 'corr') && ~isfield(m.opts, 'clust'), ...
        'opts.corr/opts.clust were not removed after merging into opts.sifted.');
    assert(isfield(m.opts, 'sifted') && m.opts.sifted.rho == 0.35, ...
        'opts.corr.threshold did not map to opts.sifted.rho.');
    assert(m.opts.sifted.flag == false, ...
        'opts.clust.flag=false did not set opts.sifted.flag=false.');

    m = ISAmigrateModel(struct('opts', struct('perf', struct('MaxMin', true))));
    assert(isfield(m.opts.perf, 'MaxPerf') && m.opts.perf.MaxPerf == true && ...
        ~isfield(m.opts.perf, 'MaxMin'), 'opts.perf.MaxMin -> opts.perf.MaxPerf rename failed.');

    m = ISAmigrateModel(struct('data', struct('bestPerformace', [1;2;3])));
    assert(isfield(m.data, 'Ybest') && isequal(m.data.Ybest, [1;2;3]) && ...
        ~isfield(m.data, 'bestPerformace'), 'model.data.bestPerformace -> model.data.Ybest rename failed.');

    m = ISAmigrateModel(struct('prelim', struct(), 'sifted', struct(), 'pilot', struct()));
    assert(isequal(m.completedStages, {'prelim','sifted','pilot'}), ...
        'completedStages was not correctly inferred from present model sub-structs.');

    lastwarn('');
    ISAmigrateModel(struct('pilot', struct('A', 1)));
    [~, warnId] = lastwarn();
    assert(strcmp(warnId, 'ISA:ISAmigrateModel:incompletePilot'), ...
        'model.pilot.A without B/C should raise an ISA:ISAmigrateModel:incompletePilot warning.');

    % LIBSVM retraining and legacy TRACE recompute need real trained
    % Z/Y/Ybin/etc. Reuse the already-built ''obj'' from the class_api
    % case above when it's usable (cheaper than rebuilding); otherwise
    % build a small fresh model here so this case stays self-contained --
    % a class_api failure (or running just this block on its own)
    % shouldn't cause a second, unrelated-looking assertion failure here
    % that hides the real one.
    if exist('obj', 'var') == 1 && isa(obj, 'InstanceSpace') && isfield(obj.model, 'trace')
        baseModel = obj.model;
    else
        fprintf('[TEST] No usable ''obj'' from class_api; building a fresh model for the LIBSVM/TRACE migration checks.\n');
        migrateObj = InstanceSpace(classCaseDir, baseOpts).build();
        baseModel = migrateObj.model;
    end

    % LIBSVM struct: a plain struct (no predict() method) under the old
    % 'svm' field name, as svmtrain() would have produced. NOTE: assigned
    % directly (pythiaLegacy.svm = {...}), not via struct('svm',{...}) --
    % the latter's cell-value form builds a struct ARRAY from the cell's
    % elements instead of a scalar struct whose field holds the cell array.
    pythiaLegacy = struct();
    pythiaLegacy.svm = repmat({struct('nr_class', 2)}, 1, numel(baseModel.data.algolabels));
    libsvmModel = baseModel;
    libsvmModel.pythia = pythiaLegacy;
    m = ISAmigrateModel(libsvmModel);
    assert(isfield(m.pythia, 'classifiers') && ~isfield(m.pythia, 'svm'), ...
        'LIBSVM-format pythia.svm was not retrained/migrated into pythia.classifiers.');
    assert(numel(m.pythia.classifiers) == numel(baseModel.data.algolabels), ...
        'Retrained pythia.classifiers has the wrong number of algorithms.');

    % #42 regression: retrainLibsvmPythia previously passed model.data.Y
    % (normalized) instead of model.data.Yraw to PYTHIA. Avg/Std_Perf_all_instances
    % are pure functions of the Y passed in (independent of which classifier
    % was trained), so a normalized-scale bug here would show up as a scale
    % mismatch against baseModel's own summary, built through the normal
    % path on the same data.
    avgPerfRetrained = cell2mat(m.pythia.summary(2:end-2, 2));
    avgPerfOriginal  = cell2mat(baseModel.pythia.summary(2:end-2, 2));
    assert(isequal(avgPerfRetrained, avgPerfOriginal), ...
        ['Retrained pythia.summary''s Avg_Perf_all_instances does not match a freshly-built ' ...
         'model''s summary on the same data -- retrainLibsvmPythia is training on the wrong ' ...
         '(normalized) scale.']);

    % #29 regression: the LIBSVM MEX-files (svmpredict/svmtrain) were
    % removed from this repository as precompiled binaries with no build
    % source in the tree. PYTHIA's eval mode must now give a clear,
    % actionable ISA:PYTHIA:noLibsvm error instead of MATLAB's generic
    % "undefined function 'svmpredict'" if it ever needs to dispatch to a
    % legacy (un-retrained) LIBSVM-format classifier struct. Uses
    % pythiaLegacy directly (not m, which was already retrained above) so
    % the isstruct(clf) legacy-LIBSVM branch is actually exercised; assumes
    % no LIBSVM MEX-files are on this MATLAB path, which holds for this
    % repository post-#29 unless one is installed separately.
    noLibsvmRaised = false;
    try
        PYTHIA(baseModel.pilot.Z, baseModel.data.Yraw, baseModel.data.Ybin, ...
            baseModel.data.Ybest, baseModel.data.algolabels, baseModel.opts.pythia, pythiaLegacy);
    catch pyErr
        noLibsvmRaised = strcmp(pyErr.identifier, 'ISA:PYTHIA:noLibsvm');
    end
    assert(noLibsvmRaised, ...
        ['PYTHIA eval mode on a legacy LIBSVM struct should raise ISA:PYTHIA:noLibsvm ' ...
         'when svmpredict is unavailable, not MATLAB''s generic undefined-function error.']);

    % Legacy DBSCAN+polyshape triangulation format: .space.area, no .measure.
    legacyTraceModel = baseModel;
    legacyTraceModel.trace = struct('space', struct('area', 1, 'polygon', []));
    m = ISAmigrateModel(legacyTraceModel);
    assert(isfield(m.trace.space, 'measure') && ~isfield(m.trace.space, 'area'), ...
        'Legacy triangulation-format model.trace was not recomputed via TRACE3.');
    assert(isfield(m.trace, 'good') && numel(m.trace.good) == numel(baseModel.data.algolabels), ...
        'Recomputed model.trace is missing per-algorithm footprints.');

    results(end).passed  = true;
    results(end).message = 'OK';
    fprintf('[TEST] Case ''isamigratemodel'' PASSED.\n');
catch ME
    results(end).passed  = false;
    results(end).message = ME.message;
    fprintf('[TEST] Case ''isamigratemodel'' FAILED: %s\n', ME.message);
end

% ---- Seed reproducibility (#41: opts.general.seed threading) ----------
% PILOT/SIFTED previously always reseeded via rng('default'), ignoring
% opts.general.seed entirely -- a user rerunning with a different seed for
% a replication/variance study got identical PILOT/SIFTED results anyway.
% Confirms a different seed changes PILOT's BFGS-multi-start result, and
% that the same seed still reproduces bit-identical results.
fprintf('\n[TEST] === Seed reproducibility: opts.general.seed threads through PILOT/SIFTED ===\n');
results(end+1).name = 'seed_reproducibility';
try
    seedOptsA = baseOpts;
    seedOptsA.general.seed = 1;
    objSeedA = InstanceSpace(classCaseDir, seedOptsA).build('stages', {'prelim', 'sifted', 'pilot'});

    seedOptsB = baseOpts;
    seedOptsB.general.seed = 2;
    objSeedB = InstanceSpace(classCaseDir, seedOptsB).build('stages', {'prelim', 'sifted', 'pilot'});

    % PILOT's numerical (BFGS) path draws its multi-start X0 from
    % rng(opts.seed,'twister'); different seeds must give different random
    % starting points and therefore (with overwhelming probability, given
    % the non-convex cost surface) a different converged A. This is the
    % robust half of the check -- SIFTED's K-means+GA feature selection can
    % legitimately converge to the same selvars across seeds on a small
    % dataset even when its own seeding is correct, so that's logged for
    % information rather than asserted.
    assert(~isequal(objSeedA.model.pilot.A, objSeedB.model.pilot.A), ...
        'Different opts.general.seed values produced identical PILOT.A -- the seed is not threading through PILOT''s BFGS multi-start.');
    if isequal(objSeedA.model.sifted.selvars, objSeedB.model.sifted.selvars)
        fprintf('[TEST] Note: SIFTED.selvars happened to match across seeds 1 and 2 (not itself a failure).\n');
    end

    % Same seed, rebuilt from scratch: must reproduce bit-identically.
    objSeedA2 = InstanceSpace(classCaseDir, seedOptsA).build('stages', {'prelim', 'sifted', 'pilot'});
    assert(isequal(objSeedA.model.pilot.A, objSeedA2.model.pilot.A) && ...
           isequal(objSeedA.model.sifted.selvars, objSeedA2.model.sifted.selvars), ...
        'Rebuilding with the same opts.general.seed changed PILOT.A/SIFTED.selvars -- reproducibility regressed.');

    results(end).passed  = true;
    results(end).message = 'OK';
    fprintf('[TEST] Case ''seed_reproducibility'' PASSED.\n');
catch ME
    results(end).passed  = false;
    results(end).message = ME.message;
    fprintf('[TEST] Case ''seed_reproducibility'' FAILED: %s\n', ME.message);
end

% ---- CLOISTER mean-centring warning (#44) ------------------------------
% CLOISTER's correlation-contradiction filter assumes mean-centred data;
% opts.auto.preproc=false or opts.norm.flag=false skips PRELIM's Box-Cox+
% Z-score step entirely, silently making the filter degenerate for any
% naturally all-positive feature. runCloister should warn in that case,
% and must not warn on the default (normalised) path.
fprintf('\n[TEST] === CLOISTER mean-centring warning: opts.norm.flag=false ===\n');
results(end+1).name = 'cloister_norm_off_warning';
try
    noNormOpts = baseOpts;
    noNormOpts.norm.flag = false;
    lastwarn('');
    InstanceSpace(classCaseDir, noNormOpts).build('stages', {'prelim', 'sifted', 'pilot', 'cloister'});
    [~, warnId] = lastwarn();
    assert(strcmp(warnId, 'ISA:InstanceSpace:cloisterNotMeanCentred'), ...
        'opts.norm.flag=false should raise ISA:InstanceSpace:cloisterNotMeanCentred before CLOISTER runs.');

    lastwarn('');
    InstanceSpace(classCaseDir, baseOpts).build('stages', {'prelim', 'sifted', 'pilot', 'cloister'});
    [~, warnId2] = lastwarn();
    assert(~strcmp(warnId2, 'ISA:InstanceSpace:cloisterNotMeanCentred'), ...
        'The default (normalised) path should not raise the cloisterNotMeanCentred warning.');

    results(end).passed  = true;
    results(end).message = 'OK';
    fprintf('[TEST] Case ''cloister_norm_off_warning'' PASSED.\n');
catch ME
    results(end).passed  = false;
    results(end).message = ME.message;
    fprintf('[TEST] Case ''cloister_norm_off_warning'' FAILED: %s\n', ME.message);
end

% ---- Per-stage required-field validation (#28) -------------------------
% checkPrereq only confirms the prerequisite STAGE completed; checkRequiredFields
% additionally confirms the specific obj.model fields that stage needs are
% actually present -- catches a hand-edited/incomplete model.mat that
% completedStages alone would wave through.
fprintf('\n[TEST] === InstanceSpace per-stage required-field validation (#28) ===\n');
results(end+1).name = 'required_fields_validation';
try
    reqObj = InstanceSpace(classCaseDir, baseOpts).build('stages', {'prelim'});
    reqObj.model.data = rmfield(reqObj.model.data, 'Ybin'); % simulate an incomplete/corrupted model
    missingFieldRaised = false;
    try
        reqObj.build('stages', {'sifted'});
    catch reqErr
        missingFieldRaised = strcmp(reqErr.identifier, 'ISA:InstanceSpace:missingField');
    end
    assert(missingFieldRaised, ...
        ['build(''stages'',{''sifted''}) with obj.model.data.Ybin removed should raise ' ...
         'ISA:InstanceSpace:missingField, not an opaque crash inside SIFTED.']);

    % Success path must be unaffected by the added check.
    reqObj2 = InstanceSpace(classCaseDir, baseOpts).build('stages', {'prelim', 'sifted'});
    assert(isequal(reqObj2.completedStages, {'prelim', 'sifted'}), ...
        'checkRequiredFields should not interfere with a normal successful build.');

    results(end).passed  = true;
    results(end).message = 'OK';
    fprintf('[TEST] Case ''required_fields_validation'' PASSED.\n');
catch ME
    results(end).passed  = false;
    results(end).message = ME.message;
    fprintf('[TEST] Case ''required_fields_validation'' FAILED: %s\n', ME.message);
end

% ---- TRACE eval-mode ngood/nbest mismatch guard (#28) -------------------
% Folded-in finding from #28: TRACE's eval mode used ngood/nbest as the
% real-vs-placeholder boundary for out.good/out.best without ever checking
% they agree, relying entirely on an unasserted cross-file convention.
fprintf('\n[TEST] === TRACE eval-mode ngood/nbest mismatch guard (#28) ===\n');
results(end+1).name = 'trace_goodbest_mismatch_guard';
try
    assert(numel(baseModel.trace.good) >= 2, ...
        'This check needs at least 2 trained algorithms to desync good/best counts.');
    mismatchedTrace = baseModel.trace;
    mismatchedTrace.good = mismatchedTrace.good(1:end-1);
    mismatchRaised = false;
    try
        TRACE(baseModel.pilot.Z, baseModel.data.Ybin, baseModel.pythia.Yhat, ...
            baseModel.data.P, baseModel.data.beta, baseModel.data.algolabels, ...
            baseModel.opts.trace, mismatchedTrace);
    catch trErr
        mismatchRaised = strcmp(trErr.identifier, 'ISA:TRACE:goodBestCountMismatch');
    end
    assert(mismatchRaised, ...
        ['TRACE eval mode with mismatched trainedTrace.good/.best counts should raise ' ...
         'ISA:TRACE:goodBestCountMismatch.']);

    results(end).passed  = true;
    results(end).message = 'OK';
    fprintf('[TEST] Case ''trace_goodbest_mismatch_guard'' PASSED.\n');
catch ME
    results(end).passed  = false;
    results(end).message = ME.message;
    fprintf('[TEST] Case ''trace_goodbest_mismatch_guard'' FAILED: %s\n', ME.message);
end

% ---- PRELIM train/eval tie-breaking consistency (#37/#38) --------------
% Before #38, InstanceSpace.evaluateTestSet reimplemented PRELIM's
% Ybin/Ybest/P/beta computation independently, using a silent
% deterministic sort() with no tie-breaking logic at all -- diverging from
% PRELIM's own random tie-break for instances with more than one
% best-performing algorithm. #38 gave PRELIM itself an eval mode
% (nargin==4) that shares this exact code with training mode, so the two
% can no longer drift apart. Verified directly here (not through a full
% build+explore round trip) with a synthetic tie, seeding rng identically
% before each call: PRELIM's tie-break is the only randomness before that
% point in the function, so an identical seed must give an identical pick
% if -- and only if -- both modes are genuinely running the same code.
fprintf('\n[TEST] === PRELIM train/eval tie-breaking consistency (#37/#38) ===\n');
results(end+1).name = 'prelim_tie_breaking_consistency';
try
    Xsynth = [1 2; 3 4; 5 6; 7 8];
    Ysynth = [1 1 5; 2 3 4; 1 5 2; 6 2 1]; % row 1: algorithms 1 and 2 tied for best (both = 1)
    prelimTestOpts = struct('MaxPerf', false, 'AbsPerf', false, 'epsilon', 0.05, ...
        'betaThreshold', 0.55, 'auto', true, 'bound', true, 'norm', true, 'iqrMultiplier', 5);

    rng(7, 'twister');
    [~, ~, trainOut] = PRELIM(Xsynth, Ysynth, prelimTestOpts);
    assert(ismember(trainOut.P(1), [1, 2]), ...
        'Training-mode tie-break for the tied instance should pick one of the tied algorithms (1 or 2).');

    rng(7, 'twister');
    [~, ~, evalOut] = PRELIM(Xsynth, Ysynth, prelimTestOpts, trainOut);
    assert(evalOut.P(1) == trainOut.P(1), ...
        ['Eval-mode PRELIM should break the same tie the same way training mode does when seeded ' ...
         'identically -- #37''s drift meant these used to disagree.']);

    results(end).passed  = true;
    results(end).message = 'OK';
    fprintf('[TEST] Case ''prelim_tie_breaking_consistency'' PASSED.\n');
catch ME
    results(end).passed  = false;
    results(end).message = ME.message;
    fprintf('[TEST] Case ''prelim_tie_breaking_consistency'' FAILED: %s\n', ME.message);
end

nCases = numel(results);

% ---- Summary ----------------------------------------------------------
fprintf('\n[TEST] ================= Summary =================\n');
nPassed = 0;
for i = 1:nCases
    if results(i).passed
        status = 'PASS';
        nPassed = nPassed + 1;
    else
        status = 'FAIL';
    end
    fprintf('[TEST] [%s] %s\n', status, results(i).name);
end
fprintf('[TEST] %d/%d cases passed.\n', nPassed, nCases);

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
