function model = ISAmigrateModel(model)
% ISAmigrateModel  Migrate a pre-v1.7 ISA model to the current field layout.
%
% Usage:
%   model = ISAmigrateModel(model)
%
% Loads a model saved by an older version of buildIS and updates field names
% to match the current PYTHIA interface:
%
%   model.pythia.svm{i}  / .knn{i}   -> model.pythia.classifiers{i}  (note: plural)
%   model.pythia.boxcosnt             -> model.pythia.param1
%   model.pythia.kscale               -> model.pythia.param2
%   opts.oracle                       -> opts.pythia  (very old models)
%
% After migration the model can be passed to PYTHIA eval mode and scriptcsv.
%
% Example:
%   m = load('model.mat');
%   m = ISAmigrateModel(m);
%   save('model.mat', '-struct', 'm');

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
