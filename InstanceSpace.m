classdef InstanceSpace
% InstanceSpace  Value-class wrapper around the ISA pipeline.
%
%   obj = InstanceSpace(rootdir)
%   obj = InstanceSpace(rootdir, opts)
%
%   Wraps PRELIM/SIFTED/PILOT/CLOISTER/PYTHIA/TRACE with stage-level
%   execution control, in-between option changes, and save()/load()
%   persistence. Because InstanceSpace is a value class (not a handle
%   class), every method that changes state returns the updated object;
%   assign the result back, e.g. obj = obj.build().
%
%   Basic usage:
%     obj = InstanceSpace(rootdir);
%     obj = obj.build();                          % run every stage
%     obj = obj.explore(testRootDir);              % evaluate on new data
%     results = obj.getResults();                  % training results
%     testResults = obj.getResults(1);              % first explore() result
%
%   Staged usage, with option changes between stages:
%     obj = InstanceSpace(rootdir);
%     obj.opts.pilot.dims = 3;
%     obj = obj.build('stages', {'prelim','sifted','pilot'});
%     figure; scatter(obj.model.pilot.Z(:,1), obj.model.pilot.Z(:,2));
%     obj.opts.pilot.alpha = 2.0;
%     obj = obj.build('stages', {'pilot'});         % re-run PILOT only
%     obj.opts.pythia.tuning = 'bayes';
%     obj = obj.build('stages', {'cloister','pythia','trace'});
%
%   Persistence:
%     obj.save();                    % writes rootdir/model.mat (-v7.3)
%     obj = InstanceSpace.load(rootdir);   % reads it back
%
%   buildIS.m and exploreIS.m are thin backward-compatibility wrappers
%   around this class; new code should prefer the class directly.

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

    properties (Access = public)
        rootdir         (1,:) char
        opts            (1,1) struct
        model           (1,1) struct
        testDirs        (1,:) cell
        testResults     (1,:) cell
        completedStages (1,:) cell
    end

    properties (Access = private, Constant)
        % Canonical execution order; build('stages',...) always runs
        % requested stages in this order regardless of how they were
        % listed, so prerequisites are satisfied within a single call.
        StageOrder = {'prelim','sifted','pilot','cloister','pythia','trace'};
        % Immediate prerequisite each stage needs already completed
        % (either from an earlier stage in the same build() call, or a
        % previous build() call on this object).
        StagePrereq = struct('prelim',   {{}}, ...
                              'sifted',   {{'prelim'}}, ...
                              'pilot',    {{'sifted'}}, ...
                              'cloister', {{'pilot'}}, ...
                              'pythia',   {{'pilot'}}, ...
                              'trace',    {{'pythia'}});
        % obj.model field name for each stage's output. Identity for
        % every stage except 'cloister', whose output runCloister()
        % stores under the shorter obj.model.cloist (matching the
        % pre-refactor buildIS.m's field name).
        StageModelField = struct('prelim',   'prelim', ...
                                  'sifted',   'sifted', ...
                                  'pilot',    'pilot', ...
                                  'cloister', 'cloist', ...
                                  'pythia',   'pythia', ...
                                  'trace',    'trace');
        % obj.model fields (dotted path, checked via hasNestedField) each
        % stage's run* method dereferences before doing any real work.
        % checkPrereq only verifies the prerequisite STAGE completed; it
        % says nothing about whether that stage's specific output fields
        % are actually present and non-empty on THIS object -- a gap that
        % matters for a model reconstructed via InstanceSpace.load(), a
        % legacy migration (ISAmigrateModel), or a hand-edited model.mat,
        % none of which re-run the stages themselves. Deliberately NOT a
        % full dependency-injection system, just the fields every
        % successful run of that stage unconditionally produces --
        % conditional-path fields (e.g. sifted's density-resubsetting
        % branch) are already isfield()-guarded at their own call site and
        % don't belong in a static required-fields list.
        StageRequiredFields = struct( ...
            'prelim',   {{}}, ...
            'sifted',   {{'data.X', 'data.Y', 'data.Ybin', 'data.featlabels'}}, ...
            'pilot',    {{'data.X', 'data.Y', 'data.featlabels'}}, ...
            'cloister', {{'data.X', 'pilot.A'}}, ...
            'pythia',   {{'pilot.Z', 'data.Yraw', 'data.Ybin', 'data.Ybest', 'data.algolabels'}}, ...
            'trace',    {{'pilot.Z', 'data.Ybin', 'pythia.Yhat', 'data.P', 'data.beta', 'data.algolabels'}});
    end

    methods (Access = public)
        function obj = InstanceSpace(rootdir, opts, requireData)
            % Reads metadata.csv presence and fills opts defaults; runs
            % no computation (spec §7.2). requireData (default true) is
            % an internal flag -- not part of the public two-argument
            % constructor call -- that InstanceSpace.load() sets false:
            % explore()-only usage (the legacy exploreIS.m flow) only
            % ever needed model.mat + metadata_test.csv on disk, so
            % load() must not depend on the original training
            % metadata.csv still being present in rootdir.
            InstanceSpace.ensurePathSetup();
            narginchk(1, 3);
            if nargin < 3 || isempty(requireData)
                requireData = true;
            end
            if ~(endsWith(rootdir, '/') || endsWith(rootdir, '\'))
                rootdir = [rootdir '/'];
            end
            obj.rootdir = rootdir;
            if requireData && ~isfile([rootdir 'metadata.csv'])
                error('ISA:InstanceSpace:missingData', ...
                    'metadata.csv not found in ''%s''.', rootdir);
            end
            if nargin < 2 || isempty(opts)
                optsfile = [rootdir 'options.json'];
                if isfile(optsfile)
                    opts = jsondecode(fileread(optsfile));
                else
                    opts = struct();
                end
            end
            opts = ISAvalidateOpts(opts);
            obj.opts             = ISAdefaults(opts);
            obj.model             = struct();
            % cell(1,0), not {}: {} is 0x0, and (1,:) property validation
            % should not be gambled on treating that as equivalent to 1x0.
            obj.testDirs         = cell(1,0);
            obj.testResults       = cell(1,0);
            obj.completedStages   = cell(1,0);
        end

        function obj = build(obj, varargin)
            % obj = obj.build() runs every stage.
            % obj = obj.build('stages', {'pilot', ...}) runs only the
            % named stages (in canonical order), erroring if a requested
            % stage's prerequisite hasn't already completed.
            p = inputParser;
            addParameter(p, 'stages', InstanceSpace.StageOrder, ...
                @(x) iscell(x) && all(ismember(x, InstanceSpace.StageOrder)));
            parse(p, varargin{:});
            toRun = InstanceSpace.StageOrder(ismember(InstanceSpace.StageOrder, p.Results.stages));

            startProcess = tic;
            fprintf('[BUILD] Root directory: %s\n', obj.rootdir);
            rng(obj.opts.general.seed, 'twister');
            if obj.opts.general.verbose
                fprintf('[BUILD] Listing options in use:\n');
                InstanceSpace.printOptions(obj.opts);
            end
            [mypool, poolOpenedHere] = obj.ensurePool();

            for i = 1:numel(toRun)
                stage = toRun{i};
                obj.checkPrereq(stage);
                obj.checkRequiredFields(stage);
                switch stage
                    case 'prelim',   obj = obj.runPrelim();
                    case 'sifted',   obj = obj.runSifted();
                    case 'pilot',    obj = obj.runPilot();
                    case 'cloister', obj = obj.runCloister();
                    case 'pythia',   obj = obj.runPythia();
                    case 'trace',    obj = obj.runTrace();
                end
                if ~ismember(stage, obj.completedStages)
                    obj.completedStages = [obj.completedStages, {stage}];
                end
                % Re-running an earlier stage invalidates every LATER
                % stage's output (e.g. re-running 'pilot' changes Z, so
                % any already-computed pythia/trace results were derived
                % from the OLD Z): drop them from completedStages and
                % remove their model fields so save()/explore()/output
                % writers can't silently operate on a model that mixes
                % the freshly (re-)run stage with stale downstream
                % results. Stages later in toRun that legitimately need
                % re-running are unaffected -- they just get re-added a
                % few iterations later, same as any other stage.
                obj = obj.invalidateDownstream(stage);
            end

            if poolOpenedHere
                fprintf('[BUILD] Closing parallel processing pool.\n');
                delete(mypool);
            end

            % Persist and write outputs only once the FULL pipeline has
            % completed (possibly across several staged build() calls, per
            % the interactive workflow in this class's help text) -- a
            % partial run (e.g. 'stages',{'prelim','sifted','pilot'}) has
            % no obj.model.pythia/trace yet, and scriptcsv/scriptpng
            % unconditionally expect every stage's output.
            if all(ismember(InstanceSpace.StageOrder, obj.completedStages))
                obj.model.opts = obj.opts;
                fprintf('[BUILD] Storing the raw MATLAB results for post-processing and/or debugging.\n');
                obj.save();
                if obj.opts.outputs.csv
                    scriptcsv(obj.model, obj.rootdir);
                    if obj.opts.outputs.web
                        scriptweb(obj.model, obj.rootdir);
                    end
                end
                if obj.opts.outputs.png
                    scriptpng(obj.model, obj.rootdir);
                end
            end
            fprintf('[BUILD] Completed in %.1f s.\n', toc(startProcess));
        end

        function obj = explore(obj, testRootDir)
            % Evaluates the trained model (obj.model) on new instances
            % from testRootDir/metadata_test.csv, using the opts frozen
            % at training time (obj.model.opts), not any opts changed on
            % obj since. Appends to obj.testDirs/obj.testResults.
            narginchk(2, 2);
            if ~(endsWith(testRootDir, '/') || endsWith(testRootDir, '\'))
                testRootDir = [testRootDir '/'];
            end
            if isempty(fieldnames(obj.model)) || ...
                    ~all(ismember(InstanceSpace.StageOrder, obj.completedStages))
                % Checking obj.model.pilot alone let a partially-built
                % object (e.g. only prelim/sifted/pilot run) reach
                % evaluateTestSet, which then fails deep inside on a
                % missing model.opts/model.pythia/model.trace field
                % instead of with a clear message here.
                error('ISA:InstanceSpace:notBuilt', ...
                    ['No fully trained model to explore -- call build() with every stage ' ...
                     '(or load a fully-built model) first.']);
            end
            datafile = [testRootDir 'metadata_test.csv'];
            if ~isfile(datafile)
                error('ISA:InstanceSpace:missingTestData', ...
                    'metadata_test.csv not found in ''%s''.', testRootDir);
            end

            startProcess = tic;
            fprintf('[EXPLORE] Root directory: %s\n', testRootDir);
            trainedModel = obj.model;
            out = InstanceSpace.evaluateTestSet(trainedModel, testRootDir);

            if trainedModel.opts.outputs.csv
                scriptcsv(out, testRootDir);
                if trainedModel.opts.outputs.web
                    scriptweb(out, testRootDir);
                end
            end
            if trainedModel.opts.outputs.png
                scriptpng(out, testRootDir);
            end

            obj.testDirs{end+1}   = testRootDir;
            obj.testResults{end+1} = out;
            fprintf('[EXPLORE] Completed in %.1f s.\n', toc(startProcess));
        end

        function results = getResults(obj, idx)
            % results = obj.getResults()    -> training results (obj.model)
            % results = obj.getResults(idx) -> the idx-th explore() result
            if nargin < 2
                results = obj.model;
            else
                if idx < 1 || idx > numel(obj.testResults)
                    error('ISA:InstanceSpace:badResultIndex', ...
                        'idx must be between 1 and %d (number of explore() calls so far).', ...
                        numel(obj.testResults));
                end
                results = obj.testResults{idx};
            end
        end

        function plot(obj, varargin)
            % plot(obj, viewName)
            % plot(obj, viewName, algoIdx)
            %
            % Interactive convenience wrapper around scriptfcn.m's drawing
            % helpers, plotting to the current figure instead of writing a
            % PNG. viewName is one of:
            %   'sources'   drawSources(Z, S)                  (needs model.data.S)
            %   'portfolio' drawPortfolioSelections(Z, P, algolabels, ...)
            %   'good'      drawBinaryPerformance(Z, Ybin(:,algoIdx), ...)
            %   'footprint' drawGoodBadFootprint(Z, good{algoIdx}, Ybin(:,algoIdx), ...)
            % algoIdx (1-based, into model.data.algolabels) is required
            % for 'good'/'footprint'.
            narginchk(2, 3);
            viewName = varargin{1};
            if ~isfield(obj.model, 'pilot')
                error('ISA:InstanceSpace:notBuilt', ...
                    'No trained model to plot -- call build() first.');
            end
            scriptfcn;
            Z = obj.model.pilot.Z;
            switch viewName
                case 'sources'
                    if ~isfield(obj.model.data, 'S')
                        error('ISA:InstanceSpace:noSources', ...
                            'model.data.S is not available (no ''source'' column in metadata.csv).');
                    end
                    drawSources(Z, obj.model.data.S);
                case 'portfolio'
                    drawPortfolioSelections(Z, obj.model.data.P, obj.model.data.algolabels, 'Best algorithm');
                case 'good'
                    algoIdx = InstanceSpace.requireAlgoIdx(varargin, obj.model.data.algolabels);
                    drawBinaryPerformance(Z, obj.model.data.Ybin(:,algoIdx), ...
                        strrep(obj.model.data.algolabels{algoIdx}, '_', ' '));
                case 'footprint'
                    algoIdx = InstanceSpace.requireAlgoIdx(varargin, obj.model.data.algolabels);
                    drawGoodBadFootprint(Z, obj.model.trace.good{algoIdx}, ...
                        obj.model.data.Ybin(:,algoIdx), strrep(obj.model.data.algolabels{algoIdx}, '_', ' '));
                otherwise
                    error('ISA:InstanceSpace:unknownView', ...
                        'Unknown plot view ''%s''. Valid views: sources, portfolio, good, footprint.', viewName);
            end
        end

        function save(obj)
            % Writes rootdir/model.mat (-v7.3, HDF5-compatible; spec §7.7).
            % Uses '-struct' to flatten obj.model's top-level fields (data,
            % prelim, sifted, pilot, cloist, pythia, trace, featsel, opts)
            % into individual MAT variables, matching the format written by
            % the pre-refactor buildIS.m so existing readers (and
            % InstanceSpace.load) keep working unchanged.
            if isempty(fieldnames(obj.model))
                error('ISA:InstanceSpace:notBuilt', ...
                    'Nothing to save -- call build() first.');
            end
            modelToSave = obj.model; %#ok<NASGU> referenced by name below via -struct
            save([obj.rootdir 'model.mat'], '-struct', 'modelToSave', '-v7.3');
        end
    end

    methods (Static)
        function obj = load(rootdir)
            % obj = InstanceSpace.load(rootdir) reads rootdir/model.mat,
            % migrates legacy field names (ISAmigrateModel) and fills any
            % option defaults absent from the saved model (ISAdefaults).
            InstanceSpace.ensurePathSetup(); % load() calls ISAmigrateModel/ISAdefaults
            % (both in utils/) directly, before the constructor call below.
            if ~(endsWith(rootdir, '/') || endsWith(rootdir, '\'))
                rootdir = [rootdir '/'];
            end
            modelfile = [rootdir 'model.mat'];
            if ~isfile(modelfile)
                error('ISA:InstanceSpace:missingModel', ...
                    'model.mat not found in ''%s''. Call build() and save() first.', rootdir);
            end
            % model.mat's top-level variables ARE model's fields (save()
            % writes them via '-struct', not a single wrapping variable),
            % so load() with an output argument reconstructs the struct
            % directly -- this also reads model.mat files written by the
            % pre-refactor buildIS.m.
            model = load(modelfile);
            model = ISAmigrateModel(model);
            model.opts = ISAdefaults(model.opts);

            obj = InstanceSpace(rootdir, model.opts, false);
            obj.model = model;
            % Infer completedStages from which stage sub-structs are
            % actually present, rather than assuming every stage ran:
            % save() is public, so a model saved after a partial
            % build('stages',{...}) call must not come back from load()
            % claiming later stages (e.g. pythia/trace) are done when
            % their fields don't exist -- that would let checkPrereq wave
            % through a build()/explore() call that then crashes deep
            % inside the missing stage instead of at the prereq check.
            obj.completedStages = InstanceSpace.StageOrder(cellfun(...
                @(s) isfield(model, InstanceSpace.StageModelField.(s)), InstanceSpace.StageOrder));
        end
    end

    methods (Access = private)
        function obj = invalidateDownstream(obj, stage)
            % Drops every stage that transitively DEPENDS ON `stage` (per
            % StagePrereq), not simply every stage later in StageOrder,
            % from completedStages, and removes its obj.model field --
            % so a re-run of an earlier stage can't leave stale
            % dependent-stage results (computed from the OLD upstream
            % output) looking valid. Called after every stage in
            % build()'s loop; a no-op when none of its dependents had
            % actually completed yet.
            %
            % StageOrder is a valid topological order but not a chain --
            % 'cloister' and 'pythia' both branch off 'pilot' and neither
            % depends on the other -- so re-running 'cloister' must NOT
            % invalidate 'pythia'/'trace' (they only depend on
            % 'pilot'/'pythia'), even though both appear later in
            % StageOrder. Found via breadth-first search over
            % StagePrereq instead of a StageOrder slice.
            dependents = {};
            frontier = {stage};
            while ~isempty(frontier)
                newFrontier = {};
                for i = 1:numel(InstanceSpace.StageOrder)
                    candidate = InstanceSpace.StageOrder{i};
                    if ismember(candidate, dependents) || ismember(candidate, frontier)
                        continue;
                    end
                    if any(ismember(InstanceSpace.StagePrereq.(candidate), frontier))
                        newFrontier{end+1} = candidate; %#ok<AGROW>
                    end
                end
                dependents = [dependents, newFrontier]; %#ok<AGROW>
                frontier = newFrontier;
            end
            for i = 1:numel(dependents)
                ds = dependents{i};
                obj.completedStages = obj.completedStages(~strcmp(obj.completedStages, ds));
                f = InstanceSpace.StageModelField.(ds);
                if isfield(obj.model, f)
                    obj.model = rmfield(obj.model, f);
                end
            end
        end

        function checkPrereq(obj, stage)
            needed = InstanceSpace.StagePrereq.(stage);
            for i = 1:numel(needed)
                if ~ismember(needed{i}, obj.completedStages)
                    error('ISA:InstanceSpace:missingPrereq', ...
                        'Stage ''%s'' requires ''%s'' to have been run first (obj.completedStages: %s).', ...
                        stage, needed{i}, strjoin(obj.completedStages, ', '));
                end
            end
        end

        function checkRequiredFields(obj, stage)
            % Field-level companion to checkPrereq (#28): checkPrereq only
            % confirms the prerequisite stage completed, not that the
            % specific obj.model fields this stage is about to dereference
            % are actually present. Raises before dispatch, naming both the
            % stage and the missing field, instead of letting it surface as
            % an opaque crash deep inside PRELIM/SIFTED/PILOT/etc.
            needed = InstanceSpace.StageRequiredFields.(stage);
            for i = 1:numel(needed)
                if ~InstanceSpace.hasNestedField(obj.model, needed{i})
                    error('ISA:InstanceSpace:missingField', ...
                        ['Stage ''%s'' requires obj.model.%s, which is missing or empty. ' ...
                         'This can happen with a hand-edited model.mat, an incomplete legacy ' ...
                         'migration (see ISAmigrateModel), or a model assembled outside the ' ...
                         'normal build() flow.'], stage, needed{i});
                end
            end
        end

        function [mypool, openedHere] = ensurePool(obj)
            % Opens a parallel pool if opts.general.parallel and none of
            % the right size already exists; reuses an existing one
            % otherwise so successive staged build() calls in the same
            % session don't pay pool-startup cost repeatedly. Only a pool
            % opened by this call is closed again at the end of build().
            mypool = gcp('nocreate');
            openedHere = false;
            if ~obj.opts.general.parallel
                return;
            end
            rightSize = ~isempty(mypool) && ...
                (~isnumeric(obj.opts.general.ncores) || mypool.NumWorkers == obj.opts.general.ncores);
            if rightSize
                return;
            end
            fprintf('[BUILD] Starting parallel processing pool.\n');
            if ~isempty(mypool)
                delete(mypool);
            end
            if isnumeric(obj.opts.general.ncores)
                mypool = parpool('local', obj.opts.general.ncores, 'SpmdEnabled', false);
            else
                mypool = parpool('local', 'SpmdEnabled', false);
            end
            openedHere = true;
        end

        function obj = runPrelim(obj)
            data = INIT(obj.rootdir, obj.opts);

            fprintf('[PRELIM] Calling PRELIM for data pre-processing.\n');
            prelimOpts = obj.opts.perf;
            prelimOpts.auto          = obj.opts.auto.preproc;
            prelimOpts.bound         = obj.opts.bound.flag;
            prelimOpts.norm          = obj.opts.norm.flag;
            prelimOpts.iqrMultiplier = obj.opts.prelim.iqrMultiplier;
            prelimOpts.nanThreshold  = obj.opts.prelim.nanThreshold;
            [data.X, data.Y, prelimOut] = PRELIM(data.X, data.Y, prelimOpts);
            data.Ybest        = prelimOut.Ybest;
            data.Ybin         = prelimOut.Ybin;
            data.P            = prelimOut.P;
            data.numGoodAlgos = prelimOut.numGoodAlgos;
            data.beta         = prelimOut.beta;

            idx = all(~data.Ybin, 1);
            if any(idx)
                warning('-> There are algorithms with no ''good'' instances. They are being removed to increase speed.');
                data.Yraw      = data.Yraw(:,~idx);
                data.Y         = data.Y(:,~idx);
                data.Ybin      = data.Ybin(:,~idx);
                data.algolabels = data.algolabels(~idx);
                if size(data.Y, 2) == 0
                    error('-> There are no ''good'' algorithms. Please verify the binary performance measure. STOPPING!')
                end
            end

            ninst = size(data.X, 1);
            fractional  = obj.opts.selvars.smallscaleflag && isfloat(obj.opts.selvars.smallscale);
            fileindexed = obj.opts.selvars.fileidxflag && isfile(obj.opts.selvars.fileidx);
            bydensity   = obj.opts.selvars.densityflag && ...
                          isfloat(obj.opts.selvars.mindistance) && ...
                          ischar(obj.opts.selvars.type);
            if fractional
                fprintf('[BUILD] Creating a small scale experiment for validation. Percentage of subset: %s%%\n', ...
                    num2str(round(100.*obj.opts.selvars.smallscale, 2)));
                state = rng;
                rng(obj.opts.general.seed, 'twister');
                aux = cvpartition(ninst, 'HoldOut', obj.opts.selvars.smallscale);
                rng(state);
                subsetIndex = aux.test;
            elseif fileindexed
                fprintf('[BUILD] Using a subset of the instances.\n');
                subsetIndex = false(size(data.X,1), 1);
                aux = table2array(readtable(obj.opts.selvars.fileidx));
                aux(aux > ninst) = [];
                subsetIndex(aux) = true;
            elseif bydensity
                fprintf('[BUILD] Creating a small scale experiment for validation based on density.\n');
                [subsetIndex, densityIsDissimilar, densityIsVISA, densityUnif] = ...
                    FILTER(data.X, data.Y, data.Ybin, obj.opts.selvars);
                subsetIndex = ~subsetIndex;
                fprintf('[BUILD] Percentage of instances retained: %s%%\n', ...
                    num2str(round(100.*mean(subsetIndex), 2)));
            else
                fprintf('[BUILD] Using the complete set of the instances.\n');
                subsetIndex = true(ninst, 1);
            end

            model_ = struct();
            model_.data = data;
            if fileindexed || fractional || bydensity
                if bydensity
                    model_.data_dense = data;
                end
                model_.data = ISAsubsetData(model_.data, subsetIndex);
            end
            model_.prelim = prelimOut;
            model_.prelim.bydensity = bydensity; % needed by the sifted stage's re-subsetting step
            if bydensity
                model_.prelim.unif = densityUnif; % feature-space uniformity of the retained subset
                % Kept for later diagnostics, not consumed by the pipeline itself: isDissimilar
                % flags instances FILTER's closeness check ruled out of the subset, isVISA flags
                % near-duplicate pairs it kept anyway (see FILTER.m's docstring).
                model_.prelim.isDissimilar = densityIsDissimilar;
                model_.prelim.isVISA = densityIsVISA;
            end
            model_.featsel.idx = 1:size(model_.data.X, 2);
            % Snapshot the pre-SIFTED feature name order here: runSifted
            % overwrites model.data.featlabels in place with the
            % SIFTED-selected subset, so this is the only place the full,
            % positionally-correct training order (post opts.selvars.feats
            % and nanThreshold drops) survives for evaluateTestSet to
            % validate metadata_test.csv's column order against.
            model_.featsel.labels = model_.data.featlabels;

            obj.model = model_;
        end

        function obj = runSifted(obj)
            nfeats = size(obj.model.data.X, 2);
            if obj.opts.sifted.flag
                fprintf('[SIFTED] Calling SIFTED for automated feature selection.\n');
                % Match the outer pipeline's final projection dimensionality
                % (spec §5.5) so feature-subset evaluation is consistent.
                siftedOpts = obj.opts.sifted;
                siftedOpts.dims = obj.opts.pilot.dims;
                [obj.model.data.X, obj.model.sifted] = SIFTED(obj.model.data.X, obj.model.data.Y, ...
                    obj.model.data.Ybin, obj.model.data.featlabels, siftedOpts);
                obj.model.data.featlabels = obj.model.data.featlabels(obj.model.sifted.selvars);
                obj.model.featsel.idx = obj.model.featsel.idx(obj.model.sifted.selvars);

                if isfield(obj.model.prelim, 'bydensity') && obj.model.prelim.bydensity
                    fprintf('[SIFTED] Creating a small scale experiment for validation based on density.\n');
                    % isDissimilar/isVISA kept for later diagnostics, not consumed by the
                    % pipeline itself (see FILTER.m's docstring).
                    [subsetIndex, obj.model.sifted.isDissimilar, obj.model.sifted.isVISA, obj.model.sifted.unif] = ...
                                         FILTER(obj.model.data_dense.X(:,obj.model.featsel.idx), ...
                                         obj.model.data_dense.Y, ...
                                         obj.model.data_dense.Ybin, ...
                                         obj.opts.selvars);
                    subsetIndex = ~subsetIndex;
                    obj.model.data = ISAsubsetData(obj.model.data_dense, subsetIndex, obj.model.featsel.idx);
                    fprintf('[SIFTED] Percentage of instances retained: %s%%\n', ...
                        num2str(round(100.*mean(subsetIndex), 2)));
                end
            else
                obj.model.sifted = struct('selvars', 1:nfeats);
            end
        end

        function obj = runPilot(obj)
            fprintf('[PILOT] Calling PILOT to find the optimal projection.\n');
            obj.model.pilot = PILOT(obj.model.data.X, obj.model.data.Y, obj.model.data.featlabels, obj.opts.pilot);
            if obj.opts.pilot.dims == 3
                fprintf('[PILOT] Finding the optimal 2D viewpoint(s) of the 3D projection.\n');
                obj.model.pilot.viewpoint = PILOTviewpoint(obj.model.pilot.Z, obj.model.data.Y, obj.opts.pilot);
            end
        end

        function obj = runCloister(obj)
            fprintf('[CLOISTER] Finding empirical bounds using CLOISTER.\n');
            % CLOISTER's correlation-contradiction filter (sign(Xedge(i,j)) ~=/==
            % sign(Xedge(i,k))) only means anything for mean-centred data -- which
            % PRELIM only produces when BOTH flags below are true (see its own
            % opts.auto && opts.norm gate). With either off, a naturally all-positive
            % feature (counts, sizes) makes the sign check degenerate for it, silently
            % making the filter a no-op for that feature rather than erroring.
            if ~(obj.opts.auto.preproc && obj.opts.norm.flag)
                warning('ISA:InstanceSpace:cloisterNotMeanCentred', ...
                    ['CLOISTER''s correlation-contradiction filter assumes mean-centred data ' ...
                     '(opts.auto.preproc and opts.norm.flag are both true), but at least one is ' ...
                     'false for this build. The sign-based filter may silently degrade for any ' ...
                     'feature that is naturally all one sign in its raw scale.']);
            end
            obj.model.cloist = CLOISTER(obj.model.data.X, obj.model.pilot.A, obj.opts.cloister);
        end

        function obj = runPythia(obj)
            fprintf('[PYTHIA] Summoning PYTHIA to train the prediction models.\n');
            obj.model.pythia = PYTHIA(obj.model.pilot.Z, obj.model.data.Yraw, obj.model.data.Ybin, ...
                obj.model.data.Ybest, obj.model.data.algolabels, obj.opts.pythia);
        end

        function obj = runTrace(obj)
            fprintf('[TRACE] Calling TRACE to perform the footprint analysis.\n');
            % Local copy: pythiaSkip is an internal signal for TRACE's
            % skip-mode detection, not a documented opts field -- must not
            % leak into opts.trace itself, since that gets persisted
            % verbatim into obj.model.opts/options.json.
            traceOpts = obj.opts.trace;
            traceOpts.pythiaSkip = obj.opts.pythia.skip;
            obj.model.trace = TRACE(obj.model.pilot.Z, obj.model.data.Ybin, obj.model.pythia.Yhat, ...
                obj.model.data.P, obj.model.data.beta, obj.model.data.algolabels, traceOpts);
        end
    end

    methods (Static, Access = private)
        function tf = hasNestedField(s, dottedPath)
            % Walks a 'a.b.c' dotted path through nested structs, used by
            % checkRequiredFields (#28). Missing at any level, or present
            % but empty at the leaf, both count as not-present -- an empty
            % array at a leaf that should hold real computed data (X, Z,
            % Yhat, ...) is exactly the "looks present, silently wrong"
            % failure mode this check exists to catch, unlike opts fields
            % (ISAvalidateOpts.getf) where an explicit [] can be a
            % deliberate, if invalid, user choice worth rejecting on its
            % own terms rather than treating as absent.
            parts = strsplit(dottedPath, '.');
            for i = 1:numel(parts)
                if ~(isstruct(s) && isfield(s, parts{i}))
                    tf = false;
                    return;
                end
                s = s.(parts{i});
            end
            tf = ~isempty(s);
        end

        function ensurePathSetup()
            % Adds this toolkit's core/output/utils/deprecated
            % subdirectories to the MATLAB path if they aren't already
            % there, so PRELIM/SIFTED/PILOT/etc. resolve regardless of
            % how the caller reached InstanceSpace -- covers buildIS,
            % exploreIS, and InstanceSpace/InstanceSpace.load used
            % directly, without requiring the caller to have run
            % startup.m first. Standalone calls to core/output/utils
            % functions themselves (bypassing InstanceSpace entirely)
            % still need startup.m -- see the note in that file.
            persistent done
            if isequal(done, true)
                return;
            end
            root = fileparts(mfilename('fullpath'));
            subdirs = {'core', 'output', 'utils', 'deprecated'};
            for i = 1:numel(subdirs)
                d = fullfile(root, subdirs{i});
                if isfolder(d)
                    addpath(d);
                end
            end
            done = true;
        end

        function printOptions(opts, prefix)
            % Compact one-line-per-setting options dump, replacing the raw
            % fieldnames()+disp() loop that used to print each top-level
            % opts section as MATLAB's verbose default struct display
            % (~80+ lines for a full options.json). Recurses into nested
            % option structs (e.g. opts.pilot, opts.sifted) so every leaf
            % setting is shown as "section.field = value".
            if nargin < 2, prefix = ''; end
            fields = fieldnames(opts);
            for i = 1:numel(fields)
                f = fields{i};
                v = opts.(f);
                key = [prefix f];
                if isstruct(v) && isscalar(v)
                    InstanceSpace.printOptions(v, [key '.']);
                else
                    fprintf('  %-28s %s\n', key, InstanceSpace.formatOptionValue(v));
                end
            end
        end

        function s = formatOptionValue(v)
            if ischar(v)
                s = ['''' v ''''];
            elseif isstring(v) && isscalar(v)
                s = ['''' char(v) ''''];
            elseif iscell(v) && all(cellfun(@ischar, v))
                s = ['{' strjoin(v, ', ') '}'];
            elseif islogical(v) && isscalar(v)
                s = mat2str(v);
            elseif isnumeric(v) && isempty(v)
                s = '[]';
            elseif isnumeric(v) && isscalar(v)
                s = num2str(v);
            elseif isnumeric(v)
                s = mat2str(v);
            else
                s = class(v);
            end
        end

        function out = evaluateTestSet(model, rootdir)
            % Ported from exploreIS.m, operating on an in-memory trained
            % model (no model.mat re-read: the caller already has it).
            [data, extra] = INIT(rootdir, model.opts, model);
            out = struct();
            out.data = data;

            fprintf('[EXPLORE] Calculating the binary measure of performance.\n');
            % model.opts.perf, flattened to PRELIM's own field names,
            % mirroring runPrelim's prelimOpts construction -- with one
            % addition kept from this function's pre-#38 form: a model
            % saved before opts.perf.MaxPerf existed (only the legacy
            % opts.perf.MaxMin, or neither) needs the same fallback
            % resolution here that ISAmigrateModel/ISAdefaults would have
            % applied for a normally-built model.
            prelimOpts = model.opts.perf;
            if ~isfield(prelimOpts, 'MaxPerf')
                if isfield(prelimOpts, 'MaxMin')
                    prelimOpts.MaxPerf = prelimOpts.MaxMin;
                else
                    warning(['Can not find parameter "MaxPerf" in the trained model. ' ...
                        'We are assuming that performance metric is needed to be minimized.']);
                    prelimOpts.MaxPerf = false;
                end
            end
            prelimOpts.auto  = model.opts.auto.preproc;
            prelimOpts.bound = model.opts.bound.flag;
            prelimOpts.norm  = model.opts.norm.flag;
            % PRELIM's eval mode recomputes Ybin/Ybest/P/beta on the new
            % (reconciled-width) Y with the exact same code -- including
            % tie-breaking -- training mode uses, then applies (not
            % re-fits) model.prelim's bounds/Box-Cox/Z-score parameters to
            % X and to the first extra.modelalgos columns of Y (#37/#38:
            % this used to be a second, independently-drifted
            % reimplementation here).
            [out.data.X, out.data.Y, prelimOut] = PRELIM(out.data.X, out.data.Y, prelimOpts, model.prelim);
            out.data.Ybest        = prelimOut.Ybest;
            out.data.Ybin         = prelimOut.Ybin;
            out.data.P            = prelimOut.P;
            out.data.numGoodAlgos = prelimOut.numGoodAlgos;
            out.data.beta         = prelimOut.beta;

            % Algorithms present in the test metadata but absent from the
            % trained model (extra.newalgos > 0) have no pre-fit
            % lambda/mu/sigma for PRELIM to apply -- fit-and-apply them
            % separately, same gating (opts.auto.preproc && opts.norm.flag)
            % as the rest of the normalisation PRELIM just did.
            if prelimOpts.auto && prelimOpts.norm && extra.newalgos > 0
                [~, out.data.Y(:,extra.modelalgos+1:extra.nalgos)] = InstanceSpace.autoNormalize( ...
                    ones(extra.ninst,1), out.data.Y(:,extra.modelalgos+1:extra.nalgos));
            end

            out.featsel.idx = model.featsel.idx;
            out.data.X = out.data.X(:,out.featsel.idx);
            out.data.featlabels = strrep(extra.featlabelsAll, 'feature_', '');
            out.data.featlabels = out.data.featlabels(model.featsel.idx);

            out.pilot.Z = out.data.X*model.pilot.A';

            out.pythia = PYTHIA(out.pilot.Z, out.data.Yraw, out.data.Ybin, out.data.Ybest, ...
                                out.data.algolabels, model.opts.pythia, model.pythia);
            out.trace = TRACE(out.pilot.Z, out.data.Ybin, out.pythia.Yhat, out.data.P, ...
                              out.data.beta, out.data.algolabels, model.opts.trace, model.trace);

            out.opts = model.opts;
        end

        function [X, Y] = autoNormalize(X, Y)
            % Fit-and-apply Box-Cox + Z-score normalisation, used only for
            % algorithms present in the test metadata but absent from the
            % training model (no pre-fit lambda/mu/sigma to reuse).
            nfeats = size(X, 2);
            nalgos = size(Y, 2);
            minX = min(X, [], 1, 'omitnan');
            X = bsxfun(@minus, X, minX) + 1;
            for i = 1:nfeats
                aux = X(:,i);
                idx = isnan(aux);
                aux = boxcox(aux(~idx));
                aux = zscore(aux);
                X(~idx,i) = aux;
            end

            minY = min(Y(:), [], 'omitnan');
            Y = (Y - minY) + eps;
            for i = 1:nalgos
                aux = Y(:,i);
                idx = isnan(aux);
                aux = boxcox(aux(~idx));
                aux = zscore(aux);
                Y(~idx,i) = aux;
            end
        end

        function algoIdx = requireAlgoIdx(args, algolabels)
            if numel(args) < 2
                error('ISA:InstanceSpace:missingAlgoIdx', ...
                    'This view needs an algorithm index: plot(obj, ''%s'', algoIdx).', args{1});
            end
            algoIdx = args{2};
            if ~(isnumeric(algoIdx) && isscalar(algoIdx) && algoIdx >= 1 && algoIdx <= numel(algolabels))
                error('ISA:InstanceSpace:badAlgoIdx', ...
                    'algoIdx must be an integer between 1 and %d.', numel(algolabels));
            end
        end
    end
end
