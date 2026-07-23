function modelOut = ISAmigrateModel(input, varargin)
% ISAmigrateModel  Migrate a pre-v1.7 ISA model to the current field layout.
%
% Two calling conventions, dispatched on the type of the first argument:
%
%   1) File-based (spec §6.4, the primary/recommended form): pass a rootdir
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
% Both paths apply the identical field migrations:
%
%   model.pythia.svm{i}  / .knn{i}   -> model.pythia.classifiers{i}  (note: plural)
%   model.pythia.boxcosnt             -> model.pythia.param1
%   model.pythia.kscale               -> model.pythia.param2
%   opts.oracle                       -> opts.pythia  (very old models)
%
% After migration the model can be passed to PYTHIA eval mode and scriptcsv.
%
% Examples:
%   ISAmigrateModel('/path/to/rootdir/');            % migrates model.mat on disk
%   m = load('model.mat'); m = ISAmigrateModel(m);   % migrates an in-memory struct

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

% ------------------------------------------------------------------
% Migrate very old opts layout: opts.oracle -> opts.pythia
if isfield(model, 'opts') && isfield(model.opts, 'oracle') && ...
        ~isfield(model.opts, 'pythia')
    model.opts.pythia = model.opts.oracle;
    model.opts = rmfield(model.opts, 'oracle');
    fprintf('ISAmigrateModel: renamed opts.oracle -> opts.pythia.\n');
end

% ------------------------------------------------------------------
% Migrate pythia struct fields.
if ~isfield(model, 'pythia')
    return;
end
p = model.pythia;

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
