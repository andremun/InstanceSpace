function opts = ISAvalidateOpts(opts)
% ISAvalidateOpts  Validate user-supplied opts fields before defaults are filled.
%   opts = ISAvalidateOpts(opts) checks the type/range of every field the
%   caller actually supplied and errors clearly (ISA:ISAvalidateOpts:*) on
%   the first invalid one, instead of letting a typo'd or out-of-range
%   value surface many stages later as a confusing crash deep inside
%   PRELIM/PILOT/PYTHIA/etc.
%
%   Deliberately validates only fields that ARE present: this runs before
%   ISAdefaults (spec §7.2), so most fields are still absent at this point
%   and are not this function's concern -- ISAdefaults supplies known-valid
%   defaults for anything missing. opts is returned unmodified; this
%   function only ever errors or passes through, it never rewrites values
%   (renaming/migrating legacy field names is ISAmigrateModel's job, not
%   this one's).

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
%
%   Simpson, C., Munoz, M.A., Kandanaarachchi, S. & Campello, R.J.G.B.
%   (2025). ISA3: A 3-dimensional expansion of Instance Space Analysis.
%   Machine Learning, 114, 240. https://doi.org/10.1007/s10994-025-06871-5
% -------------------------------------------------------------------------

if ~isstruct(opts)
    error('ISA:ISAvalidateOpts:notStruct', 'opts must be a struct; got a %s.', class(opts));
end

checkLogical('general.verbose',    getf(opts, 'general', 'verbose'));
checkLogical('general.parallel',   getf(opts, 'general', 'parallel'));
checkPosInt('general.seed',        getf(opts, 'general', 'seed'), true); % 0 allowed
% general.ncores is deliberately not type-checked to numeric here:
% InstanceSpace.ensurePool() branches on isnumeric(opts.general.ncores) and
% treats a non-numeric value (e.g. left as a placeholder, or intentionally
% not a number) as "let parpool choose its own default pool size" -- a
% supported configuration, not a misconfiguration. Only validate the range
% when it IS numeric.
ncores = getf(opts, 'general', 'ncores');
if ~isempty(ncores) && isnumeric(ncores)
    checkPosInt('general.ncores', ncores, false);
end

checkLogical('perf.MaxPerf',       getf(opts, 'perf', 'MaxPerf'));
checkLogical('perf.AbsPerf',       getf(opts, 'perf', 'AbsPerf'));
checkUnitRange('perf.epsilon',        getf(opts, 'perf', 'epsilon'));
checkUnitRange('perf.betaThreshold',  getf(opts, 'perf', 'betaThreshold'));

checkPositive('prelim.iqrMultiplier', getf(opts, 'prelim', 'iqrMultiplier'));
checkUnitRange('prelim.nanThreshold', getf(opts, 'prelim', 'nanThreshold'));

checkLogical('auto.preproc',       getf(opts, 'auto', 'preproc'));
checkLogical('bound.flag',         getf(opts, 'bound', 'flag'));
checkLogical('norm.flag',          getf(opts, 'norm', 'flag'));

checkLogical('selvars.smallscaleflag', getf(opts, 'selvars', 'smallscaleflag'));
checkUnitRange('selvars.smallscale',   getf(opts, 'selvars', 'smallscale'));
checkLogical('selvars.fileidxflag',    getf(opts, 'selvars', 'fileidxflag'));
checkText('selvars.fileidx',           getf(opts, 'selvars', 'fileidx'));
checkLogical('selvars.densityflag',    getf(opts, 'selvars', 'densityflag'));
checkPositive('selvars.mindistance',   getf(opts, 'selvars', 'mindistance'));
checkMember('selvars.type', getf(opts, 'selvars', 'type'), ...
    {'Ftr','Ftr&AP','Ftr&Good','Ftr&AP&Good'});
checkCellOfText('selvars.feats',       getf(opts, 'selvars', 'feats'));
checkCellOfText('selvars.algos',       getf(opts, 'selvars', 'algos'));

checkLogical('sifted.flag',        getf(opts, 'sifted', 'flag'));
checkUnitRange('sifted.rho',       getf(opts, 'sifted', 'rho'));
checkUnitRange('sifted.pval',      getf(opts, 'sifted', 'pval'));
checkPosInt('sifted.K',            getf(opts, 'sifted', 'K'), false);
checkPosInt('sifted.MaxIter',      getf(opts, 'sifted', 'MaxIter'), false);
checkPosInt('sifted.Replicates',   getf(opts, 'sifted', 'Replicates'), false);

checkMember('pilot.method', getf(opts, 'pilot', 'method'), {'standard','pls'});
checkMember('pilot.dims',   getf(opts, 'pilot', 'dims'), {2,3});
checkLogical('pilot.analytic',     getf(opts, 'pilot', 'analytic'));
checkPosInt('pilot.ntries',        getf(opts, 'pilot', 'ntries'), false);
checkPositive('pilot.alpha',       getf(opts, 'pilot', 'alpha'));
checkPositive('pilot.topoWeight',  getf(opts, 'pilot', 'topoWeight'), true); % 0 allowed
checkViewGroups(getf(opts, 'pilot', 'viewGroups'));

checkUnitRange('cloister.pval',            getf(opts, 'cloister', 'pval'));
checkUnitRange('cloister.corrThreshold',   getf(opts, 'cloister', 'corrThreshold'));
checkPosInt('cloister.maxFeatures',        getf(opts, 'cloister', 'maxFeatures'), false);

checkLogical('pythia.flag',        getf(opts, 'pythia', 'flag'));
checkMember('pythia.classifier', getf(opts, 'pythia', 'classifier'), ...
    {'knn','svm','tree','nb','linear','ensemble'});
checkMember('pythia.tuning', getf(opts, 'pythia', 'tuning'), {'sobol','bayes','none'});
checkPosInt('pythia.nTuningIter',  getf(opts, 'pythia', 'nTuningIter'), false);
checkPosInt('pythia.kFold',        getf(opts, 'pythia', 'kFold'), false);
checkLogical('pythia.skip',        getf(opts, 'pythia', 'skip'));
checkText('pythia.ensembleMethod', getf(opts, 'pythia', 'ensembleMethod'));

checkMember('trace.method', getf(opts, 'trace', 'method'), {'trace3','legacy'});
checkUnitRange('trace.PI',             getf(opts, 'trace', 'PI'));
checkPosInt('trace.minInstances',      getf(opts, 'trace', 'minInstances'), false);
checkUnitRange('trace.minAreaFrac',    getf(opts, 'trace', 'minAreaFrac'));
checkLogical('trace.contra',           getf(opts, 'trace', 'contra'));

checkLogical('outputs.csv',        getf(opts, 'outputs', 'csv'));
checkLogical('outputs.png',        getf(opts, 'outputs', 'png'));
checkLogical('outputs.fig',        getf(opts, 'outputs', 'fig'));
checkLogical('outputs.web',        getf(opts, 'outputs', 'web'));
end

% =========================================================================
% Helpers. Every checkX() is a no-op when passed [] (getf's sentinel for
% "field not present"), since only fields the caller actually supplied are
% validated here.

function v = getf(opts, sub, field)
if isfield(opts, sub) && isfield(opts.(sub), field)
    v = opts.(sub).(field);
else
    v = [];
end
end

function checkLogical(name, v)
if isempty(v), return; end
if ~(islogical(v) && isscalar(v))
    error('ISA:ISAvalidateOpts:notLogical', 'opts.%s must be a logical scalar (true/false); got %s.', ...
        name, describe(v));
end
end

function checkText(name, v)
if isempty(v), return; end
if ~((ischar(v)) || (isstring(v) && isscalar(v)))
    error('ISA:ISAvalidateOpts:notText', 'opts.%s must be a char or scalar string; got %s.', ...
        name, describe(v));
end
end

function checkCellOfText(name, v)
if isempty(v), return; end
if ~(iscell(v) && all(cellfun(@(x) ischar(x) || (isstring(x) && isscalar(x)), v)))
    error('ISA:ISAvalidateOpts:notCellOfText', ...
        'opts.%s must be a cell array of char/string names; got %s.', name, describe(v));
end
end

function checkPositive(name, v, zeroAllowed)
if isempty(v), return; end
if nargin < 3, zeroAllowed = false; end
if ~(isnumeric(v) && isscalar(v) && isreal(v) && ~isnan(v) && (v > 0 || (zeroAllowed && v == 0)))
    if zeroAllowed
        error('ISA:ISAvalidateOpts:notPositive', 'opts.%s must be a non-negative numeric scalar; got %s.', ...
            name, describe(v));
    else
        error('ISA:ISAvalidateOpts:notPositive', 'opts.%s must be a positive numeric scalar; got %s.', ...
            name, describe(v));
    end
end
end

function checkPosInt(name, v, zeroAllowed)
if isempty(v), return; end
if nargin < 3, zeroAllowed = false; end
checkPositive(name, v, zeroAllowed);
if ~isempty(v) && v ~= floor(v)
    error('ISA:ISAvalidateOpts:notInteger', 'opts.%s must be an integer; got %s.', name, describe(v));
end
end

function checkUnitRange(name, v)
if isempty(v), return; end
if ~(isnumeric(v) && isscalar(v) && isreal(v) && ~isnan(v) && v >= 0 && v <= 1)
    error('ISA:ISAvalidateOpts:notInUnitRange', ...
        'opts.%s must be a numeric scalar in [0, 1]; got %s.', name, describe(v));
end
end

function checkMember(name, v, validSet)
if isempty(v), return; end
if iscell(validSet) && ischar(validSet{1})
    ok = (ischar(v) || (isstring(v) && isscalar(v))) && any(strcmpi(char(v), validSet));
    optionsStr = strjoin(validSet, ''', ''');
else
    ok = isnumeric(v) && isscalar(v) && any(cellfun(@(x) isequal(v, x), validSet));
    optionsStr = strjoin(cellfun(@num2str, validSet, 'UniformOutput', false), ', ');
end
if ~ok
    error('ISA:ISAvalidateOpts:notMember', 'opts.%s must be one of {''%s''}; got %s.', ...
        name, optionsStr, describe(v));
end
end

function checkViewGroups(v)
if isempty(v), return; end
if ~iscell(v)
    error('ISA:ISAvalidateOpts:badViewGroups', ...
        'opts.pilot.viewGroups must be a cell array of algorithm index vectors; got %s.', describe(v));
end
for i = 1:numel(v)
    g = v{i};
    if ~(isnumeric(g) && isvector(g) && isreal(g) && all(g > 0) && all(g == floor(g)))
        error('ISA:ISAvalidateOpts:badViewGroups', ...
            'opts.pilot.viewGroups{%d} must be a vector of positive integer algorithm indices; got %s.', ...
            i, describe(g));
    end
end
end

function s = describe(v)
if isnumeric(v) || islogical(v)
    s = mat2str(v);
elseif ischar(v)
    s = ['''' v ''''];
elseif isstring(v) && isscalar(v)
    s = ['''' char(v) ''''];
else
    s = class(v);
end
end
