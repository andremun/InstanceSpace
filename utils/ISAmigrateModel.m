function modelOut = ISAmigrateModel(input, varargin)
% ISAmigrateModel  Migrate a legacy ISA model to the current field layout.
%
% Two calling conventions, dispatched on the type of the first argument:
%
%   1) File-based (the primary/recommended form): pass a rootdir
%      containing model.mat. The original file is backed up alongside it
%      (default name: model_legacy.mat) and the migrated model is written
%      back to model.mat in the same directory.
%
%        ISAmigrateModel(rootdir)
%        ISAmigrateModel(rootdir, 'backupSuffix', '_v1')   % -> model_v1.mat
%
%   2) In-memory: pass an already-loaded model struct and use the returned,
%      migrated struct directly; no file I/O occurs. This form exists
%      because exploreIS.m already holds the loaded model in memory (it
%      calls load(modelfile) itself to migrate-then-fill-defaults in one
%      pass) and has no reason to write the migrated model back to disk
%      and re-read it on every explore() call — the file-based form is for
%      one-time offline migration of a model.mat produced by an older
%      toolkit version, not for use on every read.
%
%        model = ISAmigrateModel(model);
%
% Both paths apply the complete legacy migration table below:
%
%   opts struct renames  opts.oracle/opts.pbldr/opts.sbound/opts.footprint
%                         -> opts.pythia/opts.pilot/opts.cloister/opts.trace
%   opts merges           opts.corr.flag/.threshold and opts.clust.flag
%                         merged into opts.sifted
%   opts.perf.MaxMin      -> opts.perf.MaxPerf
%   model.data.bestPerformace (typo) -> model.data.Ybest
%   model.pilot.A without B/C -> warning only (not expected; not auto-fixable)
%   model.pythia.svm{i} / .knn{i}   -> model.pythia.classifiers{i}
%   model.pythia.boxcosnt / .kscale -> model.pythia.param1 / .param2
%   LIBSVM struct in model.pythia   -> retrained via the current classifier
%                                      registry (opts.pythia.classifier,
%                                      default 'knn'); a LIBSVM struct has no
%                                      predict() method, so unlike the plain
%                                      field renames above it cannot simply
%                                      be relabelled
%   model.trace in the pre-refactor DBSCAN+polyshape triangulation format
%                         -> recomputed fresh via TRACE3, using
%                            model.pythia.Yhat when available (else
%                            model.data.Ybin)
%   missing model.completedStages -> inferred from which sub-structs are
%                         present (isfield only, same approach as
%                         InstanceSpace.load())
%
% After migration the model can be passed to PYTHIA eval mode and scriptcsv.
%
% Examples:
%   ISAmigrateModel('/path/to/rootdir/');            % migrates model.mat on disk
%   m = load('model.mat'); m = ISAmigrateModel(m);   % migrates an in-memory struct

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

if ischar(input) || (isstring(input) && isscalar(input))
    modelOut = migrateFromRootdir(char(input), varargin{:});
else
    if ~isempty(varargin)
        warning('ISA:ISAmigrateModel:ignoredArgs', ...
            ['Name-value arguments (e.g. ''backupSuffix'') only apply to the ' ...
             'rootdir calling form and are ignored for in-memory structs.']);
    end
    modelOut = migrateModelStruct(input);
end
end

% =========================================================================
function model = migrateFromRootdir(rootdir, varargin)
p = inputParser;
% Must be a non-empty scalar string/char: an empty suffix would make
% backupfile identical to modelfile below, so copyfile(modelfile,backupfile)
% would silently no-op (or error) instead of actually preserving a backup.
isNonEmptyText = @(x) (ischar(x) && ~isempty(x)) || ...
                      (isstring(x) && isscalar(x) && strlength(x) > 0);
addParameter(p, 'backupSuffix', '_legacy', isNonEmptyText);
parse(p, varargin{:});
backupSuffix = char(p.Results.backupSuffix);

modelfile = fullfile(rootdir, 'model.mat');
if ~isfile(modelfile)
    error('ISA:ISAmigrateModel:noModelFile', ...
        'No model.mat found in ''%s''.', rootdir);
end

model = load(modelfile);
model = migrateModelStruct(model);

[dir, base, ext] = fileparts(modelfile);
backupfile = fullfile(dir, [base backupSuffix ext]);
if strcmp(backupfile, modelfile)
    error('ISA:ISAmigrateModel:badBackupSuffix', ...
        '''backupSuffix'' (''%s'') must produce a backup filename different from model.mat.', ...
        backupSuffix);
end
if isfile(backupfile)
    % Re-running migration with the same backupSuffix would otherwise
    % silently overwrite an existing backup -- permanently losing whatever
    % pre-migration model.mat it was protecting -- instead of failing loudly.
    error('ISA:ISAmigrateModel:backupExists', ...
        ['Backup file ''%s'' already exists; refusing to overwrite it. ' ...
         'Choose a different ''backupSuffix'', or remove/rename the ' ...
         'existing backup first.'], backupfile);
end
copyfile(modelfile, backupfile);
save(modelfile, '-struct', 'model');
fprintf('ISAmigrateModel: migrated model written to ''%s''; original backed up to ''%s''.\n', ...
    modelfile, backupfile);
end

% =========================================================================
function model = migrateModelStruct(model)
if ~isstruct(model)
    return;
end

% Each migration below is independent of the others (guarded on its own
% isfield checks), so -- unlike an early-return chain -- a model missing
% one sub-struct (e.g. no model.pythia) still gets every other applicable
% migration instead of silently skipping the rest of the table.
model = migrateOptsRenames(model);
model = migrateOptsMerges(model);
model = migrateDataFieldNames(model);
model = migratePilotFields(model);
model = migratePythiaFields(model);
model = migrateTraceFields(model);
model = inferCompletedStages(model);
end

% =========================================================================
function model = migrateOptsRenames(model)
% Very old top-level opts struct renames, predating the current opts schema.
if ~isfield(model, 'opts')
    return;
end
renamePairs = {
    'oracle',    'pythia'
    'pbldr',     'pilot'
    'sbound',    'cloister'
    'footprint', 'trace'
};
for i = 1:size(renamePairs, 1)
    oldName = renamePairs{i,1};
    newName = renamePairs{i,2};
    if isfield(model.opts, oldName) && ~isfield(model.opts, newName)
        model.opts.(newName) = model.opts.(oldName);
        model.opts = rmfield(model.opts, oldName);
        fprintf('ISAmigrateModel: renamed opts.%s -> opts.%s.\n', oldName, newName);
    end
end
end

% =========================================================================
function model = migrateOptsMerges(model)
if ~isfield(model, 'opts')
    return;
end
% opts.corr.flag / opts.corr.threshold -> merged into opts.sifted. The
% correlation-selection threshold in the current schema is opts.sifted.rho
% (SIFTED.m's selection logic reads opts.rho; see ISAdefaults.m). Spec
% Section 6.4 refers to the migration target as "opts.sifted.corrThreshold",
% but no field by that literal name exists in this codebase -- rho is what
% it actually maps to.
if isfield(model.opts, 'corr')
    if ~isfield(model.opts, 'sifted')
        model.opts.sifted = struct();
    end
    if isfield(model.opts.corr, 'threshold') && ~isfield(model.opts.sifted, 'rho')
        model.opts.sifted.rho = model.opts.corr.threshold;
    end
    if isfield(model.opts.corr, 'flag') && ~model.opts.corr.flag && ...
            ~isfield(model.opts.sifted, 'flag')
        % opts.corr.flag=false meant "skip correlation-based selection",
        % which has no standalone equivalent now that correlation
        % filtering is one internal step of SIFTED rather than a separate
        % stage; disabling SIFTED entirely is the closest approximation.
        model.opts.sifted.flag = false;
    end
    model.opts = rmfield(model.opts, 'corr');
    fprintf('ISAmigrateModel: merged opts.corr into opts.sifted.\n');
end
% opts.clust.flag -> merged into opts.sifted; no direct equivalent, so a
% false flag is approximated the same way as opts.corr.flag above.
if isfield(model.opts, 'clust')
    if ~isfield(model.opts, 'sifted')
        model.opts.sifted = struct();
    end
    if isfield(model.opts.clust, 'flag') && ~model.opts.clust.flag && ...
            ~isfield(model.opts.sifted, 'flag')
        model.opts.sifted.flag = false;
    end
    model.opts = rmfield(model.opts, 'clust');
    fprintf('ISAmigrateModel: merged opts.clust into opts.sifted.\n');
end
% opts.perf.MaxMin -> opts.perf.MaxPerf
if isfield(model.opts, 'perf') && isfield(model.opts.perf, 'MaxMin') && ...
        ~isfield(model.opts.perf, 'MaxPerf')
    model.opts.perf.MaxPerf = model.opts.perf.MaxMin;
    model.opts.perf = rmfield(model.opts.perf, 'MaxMin');
    fprintf('ISAmigrateModel: renamed opts.perf.MaxMin -> opts.perf.MaxPerf.\n');
end
end

% =========================================================================
function model = migrateDataFieldNames(model)
if isfield(model, 'data') && isfield(model.data, 'bestPerformace') && ...
        ~isfield(model.data, 'Ybest')
    model.data.Ybest = model.data.bestPerformace;
    model.data = rmfield(model.data, 'bestPerformace');
    fprintf('ISAmigrateModel: renamed model.data.bestPerformace -> model.data.Ybest.\n');
end
end

% =========================================================================
function model = migratePilotFields(model)
% model.pilot.A without B/C is not expected in any production model -- B
% and C are always assigned in the same code block as A -- so there is no
% automatic fix, only a warning that the model may be corrupted or from an
% unsupported version (spec §6.4).
if isfield(model, 'pilot') && isfield(model.pilot, 'A') && ...
        (~isfield(model.pilot, 'B') || ~isfield(model.pilot, 'C'))
    warning('ISA:ISAmigrateModel:incompletePilot', ...
        ['model.pilot.A is present without B and/or C. PILOT always assigns all three ' ...
         'together, so this is unexpected; the model may be corrupted or from an ' ...
         'unsupported version. No automatic fix applied.']);
end
end

% =========================================================================
function model = migratePythiaFields(model)
if ~isfield(model, 'pythia')
    return;
end
p = model.pythia;

% LIBSVM-era classifiers are plain structs returned by svmtrain() (no
% predict() method), fundamentally incompatible with the new registry's
% ClassificationSVM objects -- unlike the old-field-name-but-already-native
% case below, these cannot simply be relabelled and must be retrained.
if isfield(p, 'svm') && ~isfield(p, 'classifiers') && ...
        iscell(p.svm) && ~isempty(p.svm) && isstruct(p.svm{1})
    model = retrainLibsvmPythia(model);
    p = model.pythia;
end

% Rename legacy classifier cell array to 'classifiers' (plural).
if isfield(p, 'svm') && ~isfield(p, 'classifiers')
    p.classifiers    = p.svm;
    p.classifierType = 'svm';
    p = rmfield(p, 'svm');
    fprintf('ISAmigrateModel: renamed pythia.svm -> pythia.classifiers (type=svm).\n');
elseif isfield(p, 'knn') && ~isfield(p, 'classifiers')
    p.classifiers    = p.knn;
    p.classifierType = 'knn';
    p = rmfield(p, 'knn');
    fprintf('ISAmigrateModel: renamed pythia.knn -> pythia.classifiers (type=knn).\n');
elseif isfield(p, 'classifier') && ~isfield(p, 'classifiers') && iscell(p.classifier)
    % Phase 4 early naming: classifier cell array stored under singular name -> plural.
    % Guard on iscell to avoid renaming the v1.7 opts string 'knn'/'svm'/etc.
    p.classifiers = p.classifier;
    p = rmfield(p, 'classifier');
    fprintf('ISAmigrateModel: renamed pythia.classifier -> pythia.classifiers.\n');
end

% Rename legacy hyperparameter fields.
if isfield(p, 'boxcosnt') && ~isfield(p, 'param1')
    p.param1 = p.boxcosnt;
    p = rmfield(p, 'boxcosnt');
    fprintf('ISAmigrateModel: renamed pythia.boxcosnt -> pythia.param1.\n');
end
if isfield(p, 'kscale') && ~isfield(p, 'param2')
    p.param2 = p.kscale;
    p = rmfield(p, 'kscale');
    fprintf('ISAmigrateModel: renamed pythia.kscale -> pythia.param2.\n');
end

% Ensure mu/sigma are present (pre-v1.7 KNN models did not z-score Z).
if ~isfield(p, 'mu') || ~isfield(p, 'sigma')
    warning('ISA:ISAmigrateModel:noZscore', ...
        ['This model has no z-score parameters (mu/sigma). ' ...
         'Test projections will not be normalised. ' ...
         'Re-train with the current PYTHIA to fix this.']);
end

model.pythia = p;
end

% =========================================================================
function model = retrainLibsvmPythia(model)
if ~isfield(model, 'pilot') || ~isfield(model.pilot, 'Z') || ...
        ~isfield(model, 'data') || ~all(isfield(model.data, {'Yraw','Ybin','Ybest','algolabels'}))
    warning('ISA:ISAmigrateModel:cannotRetrainPythia', ...
        ['model.pythia holds LIBSVM-format classifiers, which have no predict() method ' ...
         'and cannot simply be renamed -- but the fields needed to retrain them ' ...
         '(model.pilot.Z and model.data.Yraw/Ybin/Ybest/algolabels) are missing. Skipping ' ...
         'retraining: field names are still renamed (e.g. pythia.svm -> pythia.classifiers), ' ...
         'but the classifiers themselves remain legacy LIBSVM structs. The LIBSVM MEX-files ' ...
         'are no longer bundled with this repository (see README.md''s LIBSVM section) -- ' ...
         'evaluating this model as-is will error unless you obtain LIBSVM yourself ' ...
         '(https://www.csie.ntu.edu.tw/~cjlin/libsvm/). Re-run build() with the original ' ...
         'training data to retrain with the current registry instead (recommended).']);
    return;
end
pyOpts = struct();
if isfield(model, 'opts') && isfield(model.opts, 'pythia')
    pyOpts = model.opts.pythia;
end
validClassifiers = {'knn','svm','tree','nb','linear','ensemble'};
% Accept char or scalar string (ISAvalidateOpts and the rest of the opts
% schema allow both for text fields), not char only -- a valid
% opts.pythia.classifier supplied as a string would otherwise be silently
% discarded and replaced with the 'knn' fallback below. isfield() is
% checked first: pyOpts.classifier would error on a struct that doesn't
% have the field at all.
hasClassifier = isfield(pyOpts, 'classifier');
isTextClassifier = hasClassifier && (ischar(pyOpts.classifier) || ...
    (isstring(pyOpts.classifier) && isscalar(pyOpts.classifier)));
if isTextClassifier && any(strcmpi(char(pyOpts.classifier), validClassifiers))
    pyOpts.classifier = lower(char(pyOpts.classifier)); % preserve valid intent regardless of case/type
else
    pyOpts.classifier = 'knn'; % spec §1.2 default for LIBSVM migration
end
full = ISAdefaults(struct('pythia', pyOpts));
pyOpts = full.pythia;
fprintf(['ISAmigrateModel: retraining LIBSVM-era pythia classifiers using the native ' ...
    '''%s'' registry entry (opts.pythia.classifier).\n'], pyOpts.classifier);
model.pythia = PYTHIA(model.pilot.Z, model.data.Yraw, model.data.Ybin, model.data.Ybest, ...
    model.data.algolabels, pyOpts);
end

% =========================================================================
function model = migrateTraceFields(model)
if ~isfield(model, 'trace') || ~isTraceLegacyFormat(model.trace)
    return;
end
if ~isfield(model, 'pilot') || ~isfield(model.pilot, 'Z') || ...
        ~isfield(model, 'data') || ~all(isfield(model.data, {'Ybin','P','beta','algolabels'}))
    warning('ISA:ISAmigrateModel:cannotRecomputeTrace', ...
        ['model.trace is in the pre-refactor DBSCAN+polyshape triangulation format, but ' ...
         'the fields needed to recompute it via TRACE3 (model.pilot.Z and model.data.' ...
         'Ybin/P/beta/algolabels) are missing. Leaving model.trace unmigrated; re-run ' ...
         'build() to produce current-format footprints.']);
    return;
end
if isfield(model, 'pythia') && isfield(model.pythia, 'Yhat')
    Yhat = model.pythia.Yhat;
else
    % Spec §6.4 fallback: no PYTHIA predictions available to filter by, so
    % TRACE (via its own pythiaAvailable check) falls back to true labels.
    Yhat = [];
end
traceOpts = struct();
if isfield(model, 'opts') && isfield(model.opts, 'trace')
    traceOpts = model.opts.trace;
end
traceOpts.method = 'trace3';
full = ISAdefaults(struct('trace', traceOpts));
traceOpts = full.trace;
fprintf(['ISAmigrateModel: recomputing model.trace footprints via TRACE3 ' ...
    '(was pre-refactor DBSCAN+polyshape triangulation format).\n']);
model.trace = TRACE(model.pilot.Z, model.data.Ybin, Yhat, model.data.P, model.data.beta, ...
    model.data.algolabels, traceOpts);
end

function tf = isTraceLegacyFormat(traceStruct)
% Old format (TRACE_legacy.m's TRACEbuild): footprint.space.area, no
% .measure. model.trace.space is always populated (built over every
% instance, never a degenerate "throw" case), unlike individual
% good{i}/best{i} entries which can be empty for algorithms with too few
% good instances -- so space is the more reliable format signal.
tf = isfield(traceStruct, 'space') && isstruct(traceStruct.space) && ...
    isfield(traceStruct.space, 'area') && ~isfield(traceStruct.space, 'measure');
end

% =========================================================================
function model = inferCompletedStages(model)
if isfield(model, 'completedStages')
    return;
end
% Mirrors InstanceSpace.load()'s own stage-name -> model-field-name mapping
% (identity for every stage except 'cloister', whose output is stored under
% the shorter model.cloist -- matching the pre-refactor buildIS.m's field
% name). Duplicated here rather than shared because ISAmigrateModel.m
% operates on plain structs, with no dependency on the InstanceSpace class.
stageOrder = {'prelim','sifted','pilot','cloister','pythia','trace'};
stageField = struct('prelim',   'prelim', ...
                     'sifted',   'sifted', ...
                     'pilot',    'pilot', ...
                     'cloister', 'cloist', ...
                     'pythia',   'pythia', ...
                     'trace',    'trace');
present = false(1, numel(stageOrder));
for i = 1:numel(stageOrder)
    present(i) = isfield(model, stageField.(stageOrder{i}));
end
model.completedStages = stageOrder(present);
fprintf('ISAmigrateModel: inferred completedStages = {%s} from present model sub-structs.\n', ...
    strjoin(model.completedStages, ', '));
end
