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
%   Simpson, C., Munoz, M.A., Kandanaarachchi, S. & Campello, R.J.G.B.
%   (2025). ISA3: A 3-dimensional expansion of Instance Space Analysis.
%   Machine Learning, 114, 240. https://doi.org/10.1007/s10994-025-06871-5
%
%   Smith-Miles, K. & Munoz, M.A. (2023). Instance Space Analysis for
%   Algorithm Testing. ACM Computing Surveys, 55(12), Article 255.
%   https://doi.org/10.1145/3572895
% -------------------------------------------------------------------------
classdef InstanceSpace
% InstanceSpace  Value-class wrapper around the ISA pipeline (spec §7).
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
%   around this class (spec §7.5); new code should prefer the class
%   directly.

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
    end

    methods (Access = public)
        function obj = InstanceSpace(rootdir, opts)
            % Reads metadata.csv presence and fills opts defaults; runs
            % no computation (spec §7.2).
            InstanceSpace.ensurePathSetup();
            narginchk(1, 2);
            if ~(endsWith(rootdir, '/') || endsWith(rootdir, '\'))
                rootdir = [rootdir '/'];
            end
            obj.rootdir = rootdir;
            if ~isfile([rootdir 'metadata.csv'])
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
            if isempty(fieldnames(obj.model)) || ~isfield(obj.model, 'pilot')
                error('ISA:InstanceSpace:notBuilt', ...
                    'No trained model to explore -- call build() (or load an existing one) first.');
            end
            datafile = [testRootDir 'metadata_test.csv'];
            if ~isfile(datafile)
                error('ISA:InstanceSpace:missingTestData', ...
                    'metadata_test.csv not found in ''%s''.', testRootDir);
            end

            startProcess = tic;
            fprintf('[EXPLORE] Root directory: %s\n', testRootDir);
            trainedModel = obj.model;
            out = InstanceSpace.evaluateTestSet(trainedModel, datafile);

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

            obj = InstanceSpace(rootdir, model.opts);
            obj.model = model;
            % Infer completedStages from which stage sub-structs are
            % actually present, rather than assuming every stage ran:
            % save() is public, so a model saved after a partial
            % build('stages',{...}) call must not come back from load()
            % claiming later stages (e.g. pythia/trace) are done when
            % their fields don't exist -- that would let checkPrereq wave
            % through a build()/explore() call that then crashes deep
            % inside the missing stage instead of at the prereq check.
            obj.completedStages = InstanceSpace.StageOrder(...
                cellfun(@(s) isfield(model, s), InstanceSpace.StageOrder));
        end
    end

    methods (Access = private)
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
            delete(mypool);
            if isnumeric(obj.opts.general.ncores)
                mypool = parpool('local', obj.opts.general.ncores, 'SpmdEnabled', false);
            else
                mypool = parpool('local', 'SpmdEnabled', false);
            end
            openedHere = true;
        end

        function obj = runPrelim(obj)
            fprintf('[BUILD] Loading metadata.csv.\n');
            datafile = [obj.rootdir 'metadata.csv'];
            Xbar = readtable(datafile);
            varlabels = Xbar.Properties.VariableNames;
            isname   = strcmpi(varlabels, 'instances');
            isfeat   = strncmpi(varlabels, 'feature_', 8);
            isalgo   = strncmpi(varlabels, 'algo_', 5);
            issource = strcmpi(varlabels, 'source');

            data = struct();
            data.instlabels = Xbar{:,isname};
            if isnumeric(data.instlabels)
                data.instlabels = num2cell(data.instlabels);
                data.instlabels = cellfun(@(x) num2str(x), data.instlabels, 'UniformOutput', false);
            end
            if any(issource)
                data.S = categorical(Xbar{:,issource});
            end
            data.X = Xbar{:,isfeat};
            data.Y = Xbar{:,isalgo};

            data.featlabels = varlabels(isfeat);
            if isfield(obj.opts, 'selvars') && isfield(obj.opts.selvars, 'feats')
                msg = 'Using the following features: ';
                isselfeat = false(1, length(data.featlabels));
                for i = 1:length(obj.opts.selvars.feats)
                    isselfeat = isselfeat | strcmp(data.featlabels, obj.opts.selvars.feats{i});
                    msg = [msg obj.opts.selvars.feats{i} ' ']; %#ok<AGROW>
                end
                fprintf('[BUILD] %s\n', msg);
                data.X = data.X(:,isselfeat);
                data.featlabels = data.featlabels(isselfeat);
            end

            data.algolabels = varlabels(isalgo);
            if isfield(obj.opts, 'selvars') && isfield(obj.opts.selvars, 'algos')
                msg = 'Using the following algorithms: ';
                isselalgo = false(1, length(data.algolabels));
                for i = 1:length(obj.opts.selvars.algos)
                    isselalgo = isselalgo | strcmp(data.algolabels, obj.opts.selvars.algos{i});
                    msg = [msg obj.opts.selvars.algos{i} ' ']; %#ok<AGROW>
                end
                fprintf('[BUILD] %s\n', msg);
                data.Y = data.Y(:,isselalgo);
                data.algolabels = data.algolabels(isselalgo);
            end

            idx = all(isnan(data.X), 2) | all(isnan(data.Y), 2);
            if any(idx)
                warning('-> There are instances with too many missing values. They are being removed to increase speed.');
                data.X = data.X(~idx,:);
                data.Y = data.Y(~idx,:);
                data.instlabels = data.instlabels(~idx);
                if isfield(data, 'S')
                    data.S = data.S(~idx);
                end
            end
            idx = mean(isnan(data.X), 1) >= obj.opts.prelim.nanThreshold;
            if any(idx)
                warning('-> There are features with too many missing values. They are being removed to increase speed.');
                data.X = data.X(:,~idx);
                data.featlabels = data.featlabels(~idx);
            end
            ninst = size(data.X, 1);
            nuinst = size(unique(data.X, 'rows'), 1);
            if nuinst/ninst < 0.5
                warning('-> There are too many repeated instances. It is unlikely that this run will produce good results.');
            end

            data.Xraw = data.X;
            data.Yraw = data.Y;
            data.featlabels = strrep(data.featlabels, 'feature_', '');
            data.algolabels = strrep(data.algolabels, 'algo_', '');

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
                rng('default');
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
                [subsetIndex, ~, ~, densityUnif] = FILTER(data.X, data.Y, data.Ybin, obj.opts.selvars);
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
            end
            model_.featsel.idx = 1:size(model_.data.X, 2);

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
                    [subsetIndex, ~, ~, obj.model.sifted.unif] = FILTER(obj.model.data_dense.X(:,obj.model.featsel.idx), ...
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

        function out = evaluateTestSet(model, datafile)
            % Ported from exploreIS.m, operating on an in-memory trained
            % model (no model.mat re-read: the caller already has it).
            fprintf('[EXPLORE] Loading metadata_test.csv.\n');
            Xbar = readtable(datafile);
            varlabels = Xbar.Properties.VariableNames;
            isname   = strcmpi(varlabels, 'instances');
            isfeat   = strncmpi(varlabels, 'feature_', 8);
            isalgo   = strncmpi(varlabels, 'algo_', 5);
            issource = strcmpi(varlabels, 'source');

            out = struct();
            out.data.instlabels = Xbar{:,isname};
            if isnumeric(out.data.instlabels)
                out.data.instlabels = num2cell(out.data.instlabels);
                out.data.instlabels = cellfun(@(x) num2str(x), out.data.instlabels, 'UniformOutput', false);
            end
            if any(issource)
                out.data.S = categorical(Xbar{:,issource});
            end
            out.data.X = Xbar{:,isfeat};
            out.data.Y = Xbar{:,isalgo};
            [ninst, nalgos] = size(out.data.Y);

            % Mirror runPrelim's opts.selvars.feats restriction (applied at
            % build time, before PRELIM computed model.prelim.hibound/
            % lobound/minX/etc, all of which are sized to the RESTRICTED
            % feature count): without this, out.data.X keeps every raw
            % feature_ column from metadata_test.csv, and the bound-clipping
            % below fails with a dimension mismatch against those
            % restricted-size arrays -- or, if it happened not to error,
            % model.featsel.idx would silently select the wrong columns
            % positionally instead of the intended named features.
            featlabelsAll = varlabels(isfeat);
            if isfield(model.opts, 'selvars') && isfield(model.opts.selvars, 'feats')
                isselfeat = false(1, length(featlabelsAll));
                for i = 1:length(model.opts.selvars.feats)
                    isselfeat = isselfeat | strcmp(featlabelsAll, model.opts.selvars.feats{i});
                end
                out.data.X = out.data.X(:,isselfeat);
                featlabelsAll = featlabelsAll(isselfeat);
            end
            % Validate against model.prelim.lambdaX (one Box-Cox lambda per
            % feature PRELIM actually fit at training time, i.e. after both
            % opts.selvars.feats AND any opts.prelim.nanThreshold-triggered
            % column drops -- the latter isn't mirrored above, since which
            % columns those were isn't retained anywhere in the model). A
            % count mismatch here means metadata_test.csv doesn't match the
            % training feature set, and would otherwise surface many rows
            % down as an opaque bsxfun dimension-mismatch error.
            if numel(featlabelsAll) ~= numel(model.prelim.lambdaX)
                error('ISA:InstanceSpace:featureCountMismatch', ...
                    ['metadata_test.csv has %d feature column(s) after applying ' ...
                     'opts.selvars.feats, but the trained model expects %d ' ...
                     '(some training features may have been dropped for exceeding ' ...
                     'opts.prelim.nanThreshold). Check that metadata_test.csv has the ' ...
                     'same feature_ columns as metadata.csv.'], ...
                    numel(featlabelsAll), numel(model.prelim.lambdaX));
            end

            % Reconcile test algorithms against the trained model's: known
            % algorithms line up by name, unseen ones are appended as new
            % columns (NaN for training-only algorithms).
            out.data.algolabels = strrep(varlabels(isalgo), 'algo_', '');
            algoexist = zeros(1, nalgos);
            for ii = 1:nalgos
                aux = find(strcmpi(strtrim(out.data.algolabels{ii}), strtrim(model.data.algolabels)));
                if ~isempty(aux)
                    algoexist(ii) = aux;
                end
            end
            newalgos  = sum(algoexist == 0);
            modelalgos = length(model.data.algolabels);
            Yaux  = NaN + ones(ninst, modelalgos+newalgos);
            lblaux = model.data.algolabels;
            acc = modelalgos + 1;
            for ii = 1:nalgos
                if algoexist(ii) == 0
                    Yaux(:,acc) = out.data.Y(:,ii);
                    lblaux(:,acc) = out.data.algolabels(ii);
                    acc = acc + 1;
                else
                    Yaux(:,algoexist(ii)) = out.data.Y(:,ii);
                end
            end
            out.data.Y = Yaux;
            out.data.algolabels = lblaux;
            nalgos = size(out.data.Y, 2);

            out.data.Xraw = out.data.X;
            out.data.Yraw = out.data.Y;

            fprintf('[EXPLORE] Calculating the binary measure of performance.\n');
            msg = 'An algorithm is good if its performance is ';
            MaxPerf = false;
            if isfield(model.opts.perf, 'MaxPerf')
                MaxPerf = model.opts.perf.MaxPerf;
            elseif isfield(model.opts.perf, 'MaxMin')
                MaxPerf = model.opts.perf.MaxMin;
            else
                warning('Can not find parameter "MaxPerf" in the trained model. We are assuming that performance metric is needed to be minimized.');
            end
            if MaxPerf
                Yaux = out.data.Y;
                Yaux(isnan(Yaux)) = -Inf;
                [rankPerf, rankAlgo] = sort(Yaux, 2, 'descend');
                out.data.Ybest = rankPerf(:,1);
                out.data.P = rankAlgo(:,1);
                if model.opts.perf.AbsPerf
                    out.data.Ybin = out.data.Y >= model.opts.perf.epsilon;
                    msg = [msg 'higher than ' num2str(model.opts.perf.epsilon)];
                else
                    out.data.Ybest(out.data.Ybest == 0) = eps;
                    out.data.Y(out.data.Y == 0) = eps;
                    out.data.Y = 1 - bsxfun(@rdivide, out.data.Y, out.data.Ybest);
                    out.data.Ybin = (1 - bsxfun(@rdivide, Yaux, out.data.Ybest)) <= model.opts.perf.epsilon;
                    msg = [msg 'within ' num2str(round(100.*model.opts.perf.epsilon)) '% of the best.'];
                end
            else
                Yaux = out.data.Y;
                Yaux(isnan(Yaux)) = Inf;
                [rankPerf, rankAlgo] = sort(Yaux, 2, 'ascend');
                out.data.Ybest = rankPerf(:,1);
                out.data.P = rankAlgo(:,1);
                if model.opts.perf.AbsPerf
                    out.data.Ybin = out.data.Y <= model.opts.perf.epsilon;
                    msg = [msg 'less than ' num2str(model.opts.perf.epsilon)];
                else
                    out.data.Ybest(out.data.Ybest == 0) = eps;
                    out.data.Y(out.data.Y == 0) = eps;
                    out.data.Y = bsxfun(@rdivide, out.data.Y, out.data.Ybest) - 1;
                    out.data.Ybin = (bsxfun(@rdivide, Yaux, out.data.Ybest) - 1) <= model.opts.perf.epsilon;
                    msg = [msg 'within ' num2str(round(100.*model.opts.perf.epsilon)) '% of the best.'];
                end
            end
            fprintf('[EXPLORE] %s\n', msg);
            out.data.numGoodAlgos = sum(out.data.Ybin, 2);
            out.data.beta = out.data.numGoodAlgos > model.opts.perf.betaThreshold*nalgos;

            if model.opts.auto.preproc && model.opts.bound.flag
                fprintf('[EXPLORE] Auto-pre-processing. Bounding outliers, scaling and normalizing the data.\n');
                fprintf('[EXPLORE] Removing extreme outliers from the feature values.\n');
                himask = bsxfun(@gt, out.data.X, model.prelim.hibound);
                lomask = bsxfun(@lt, out.data.X, model.prelim.lobound);
                % Out-of-distribution warning (spec §7.5): more than 5% of
                % test instances clipped to the training bounds suggests
                % they may fall outside the training distribution.
                fracClipped = mean(any(himask | lomask, 2));
                if fracClipped > 0.05
                    warning('ISA:InstanceSpace:outOfDistribution', ...
                        ['%.1f%% of test instances were clipped to the training feature bounds; ' ...
                         'they may be outside the training distribution. Consider retraining with ' ...
                         'a combined dataset.'], 100*fracClipped);
                end
                out.data.X = out.data.X.*~(himask | lomask) + bsxfun(@times, himask, model.prelim.hibound) + ...
                                                              bsxfun(@times, lomask, model.prelim.lobound);
            end

            if model.opts.auto.preproc && model.opts.norm.flag
                fprintf('[EXPLORE] Auto-normalizing the data.\n');
                out.data.X = bsxfun(@minus, out.data.X, model.prelim.minX) + 1;
                out.data.X(~isnan(out.data.X) & out.data.X < 1) = 1;
                for ii = 1:length(model.prelim.lambdaX)
                    lambda = model.prelim.lambdaX(ii);
                    x = out.data.X(:,ii);
                    idx = ~isnan(x);
                    if abs(lambda) < 1e-10
                        x(idx) = log(x(idx));
                    else
                        x(idx) = (x(idx).^lambda - 1) ./ lambda;
                    end
                    out.data.X(:,ii) = x;
                end
                out.data.X = bsxfun(@rdivide, bsxfun(@minus, out.data.X, model.prelim.muX), model.prelim.sigmaX);

                out.data.Y = (out.data.Y - model.prelim.minY) + eps;
                out.data.Y(out.data.Y <= 0) = eps;
                for ii = 1:modelalgos
                    lambda = model.prelim.lambdaY(ii);
                    y = out.data.Y(:,ii);
                    idx = ~isnan(y);
                    if abs(lambda) < 1e-10
                        y(idx) = log(y(idx));
                    else
                        y(idx) = (y(idx).^lambda - 1) ./ lambda;
                    end
                    out.data.Y(:,ii) = y;
                end
                out.data.Y(:,1:modelalgos) = bsxfun(@rdivide, bsxfun(@minus, out.data.Y(:,1:modelalgos), ...
                    model.prelim.muY), model.prelim.sigmaY);
                if newalgos > 0
                    [~, out.data.Y(:,modelalgos+1:nalgos)] = InstanceSpace.autoNormalize( ...
                        ones(ninst,1), out.data.Y(:,modelalgos+1:nalgos));
                end
            end
            if ~isreal(out.data.X)
                error('ISA:InstanceSpace:complexX', ...
                    'Feature matrix X is complex after normalisation. Check test data range vs training data.');
            end

            out.featsel.idx = model.featsel.idx;
            out.data.X = out.data.X(:,out.featsel.idx);
            out.data.featlabels = strrep(featlabelsAll, 'feature_', '');
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
