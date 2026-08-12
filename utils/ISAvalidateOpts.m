function opts = ISAvalidateOpts(opts)
% ISAvalidateOpts  Validate user-supplied opts fields before defaults are filled.
%   opts = ISAvalidateOpts(opts) checks the type/range of every RECOGNISED
%   opts field the caller actually supplied (the fixed set of fields this
%   function knows about; see the body) and errors clearly
%   (ISA:ISAvalidateOpts:*) on the first invalid one, instead of letting
%   an out-of-range value surface many stages later as a confusing crash
%   deep inside PRELIM/PILOT/PYTHIA/etc. An unrecognised field name (e.g.
%   a typo like opts.piyhia.classifier) is NOT flagged -- it passes
%   through silently, exactly like an unset one, since this function has
%   no way to distinguish "not a real option" from "a future option it
%   doesn't know about yet".
%
%   Deliberately validates only fields that ARE present: this runs before
%   ISAdefaults, so most fields are still absent at this point
%   and are not this function's concern -- ISAdefaults supplies known-valid
%   defaults for anything missing. "Present" is tracked explicitly (getf()
%   returns a presence flag alongside the value), not inferred from
%   isempty(v): a field explicitly supplied as opts.general.parallel = []
%   IS present and must be rejected as invalid, not silently skipped as if
%   absent -- ISAdefaults only checks isfield(), so it would never replace
%   that [] with the proper default, and later code expecting a logical
%   scalar would fail far from the actual mistake. opts is returned
%   unmodified; this function only ever errors or passes through, it
%   never rewrites values (renaming/migrating legacy field names is
%   ISAmigrateModel's job, not this one's).

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

if ~isstruct(opts)
    error('ISA:ISAvalidateOpts:notStruct', 'opts must be a struct; got a %s.', class(opts));
end

checkLogical(opts, 'general', 'verbose');
checkLogical(opts, 'general', 'parallel');
checkPosInt(opts, 'general', 'seed', true); % 0 allowed
% general.ncores is deliberately not type-checked to numeric via the
% checkPosInt wrapper: InstanceSpace.ensurePool() branches on
% isnumeric(opts.general.ncores) and treats a non-numeric value as "let
% parpool choose its own default pool size" -- a supported configuration,
% not a misconfiguration. Only validate the range when it IS numeric.
[ncores, ncoresPresent] = getf(opts, 'general', 'ncores');
if ncoresPresent && isnumeric(ncores)
    checkPosIntValue('general.ncores', ncores, false);
end

checkLogical(opts, 'perf', 'MaxPerf');
checkLogical(opts, 'perf', 'AbsPerf');
% epsilon is only a [0,1] fraction when performance is relative (compared
% against a normalized ratio to the best algorithm); under AbsPerf=true it
% is compared directly against the raw performance measure (PRELIM.m), so
% it carries that measure's units and can be any real number.
[absPerf, absPerfPresent] = getf(opts, 'perf', 'AbsPerf');
if ~(absPerfPresent && absPerf)
    checkUnitRange(opts, 'perf', 'epsilon');
end
checkUnitRange(opts, 'perf', 'betaThreshold');

checkPositive(opts, 'prelim', 'iqrMultiplier');
checkUnitRange(opts, 'prelim', 'nanThreshold');

checkLogical(opts, 'auto', 'preproc');
checkLogical(opts, 'bound', 'flag');
checkLogical(opts, 'norm', 'flag');

checkLogical(opts, 'selvars', 'smallscaleflag');
checkUnitRange(opts, 'selvars', 'smallscale');
checkLogical(opts, 'selvars', 'fileidxflag');
checkText(opts, 'selvars', 'fileidx');
checkLogical(opts, 'selvars', 'densityflag');
checkPositive(opts, 'selvars', 'mindistance');
checkMember(opts, 'selvars', 'type', {'Ftr','Ftr&AP','Ftr&Good','Ftr&AP&Good'});
checkCellOfText(opts, 'selvars', 'feats');
checkCellOfText(opts, 'selvars', 'algos');

checkLogical(opts, 'sifted', 'flag');
checkUnitRange(opts, 'sifted', 'rho');
checkUnitRange(opts, 'sifted', 'pval');
checkPosInt(opts, 'sifted', 'K', false);
checkPosInt(opts, 'sifted', 'MaxIter', false);
checkPosInt(opts, 'sifted', 'Replicates', false);

checkMember(opts, 'pilot', 'method', {'standard','pls'});
checkMember(opts, 'pilot', 'dims', {2,3});
checkLogical(opts, 'pilot', 'analytic');
checkPosInt(opts, 'pilot', 'ntries', false);
checkPositive(opts, 'pilot', 'alpha');
checkPositive(opts, 'pilot', 'topoWeight', true); % 0 allowed
checkViewGroups(opts, 'pilot', 'viewGroups');

checkUnitRange(opts, 'cloister', 'pval');
checkUnitRange(opts, 'cloister', 'corrThreshold');
checkPosInt(opts, 'cloister', 'maxFeatures', false);

checkLogical(opts, 'pythia', 'flag');
checkMember(opts, 'pythia', 'classifier', {'knn','svm','tree','nb','linear','ensemble'});
checkMember(opts, 'pythia', 'tuning', {'sobol','bayes','none'});
checkPosInt(opts, 'pythia', 'nTuningIter', false);
checkPosInt(opts, 'pythia', 'kFold', false);
checkLogical(opts, 'pythia', 'skip');
checkText(opts, 'pythia', 'ensembleMethod');

checkMember(opts, 'trace', 'method', {'trace3','legacy'});
checkUnitRange(opts, 'trace', 'PI');
checkPosInt(opts, 'trace', 'minInstances', false);
checkUnitRange(opts, 'trace', 'minAreaFrac');
checkLogical(opts, 'trace', 'contra');

checkLogical(opts, 'outputs', 'csv');
checkLogical(opts, 'outputs', 'png');
checkLogical(opts, 'outputs', 'fig');
checkLogical(opts, 'outputs', 'web');
end

% =========================================================================
% Helpers.
%
% getf() returns [value, present]: present distinguishes "field absent"
% from "field present with an empty/invalid value" -- the two are NOT the
% same thing, so every checkX() wrapper below gates on present, never on
% isempty(v). The wrappers do the opts.(sub).(field) lookup and presence
% gating; the matching checkXValue() functions do the actual type/range
% check on an already-known-present value and are reused directly (e.g.
% general.ncores above) wherever a wrapper's all-or-nothing gating isn't
% quite right.

function [v, present] = getf(opts, sub, field)
if ~isfield(opts, sub)
    v = [];
    present = false;
    return;
end
if ~isstruct(opts.(sub))
    error('ISA:ISAvalidateOpts:notStruct', 'opts.%s must be a struct; got %s.', sub, describe(opts.(sub)));
end
if isfield(opts.(sub), field)
    v = opts.(sub).(field);
    present = true;
else
    v = [];
    present = false;
end
end

function checkLogical(opts, sub, field)
[v, present] = getf(opts, sub, field);
if ~present, return; end
checkLogicalValue([sub '.' field], v);
end

function checkLogicalValue(name, v)
if ~(islogical(v) && isscalar(v))
    error('ISA:ISAvalidateOpts:notLogical', 'opts.%s must be a logical scalar (true/false); got %s.', ...
        name, describe(v));
end
end

function checkText(opts, sub, field)
[v, present] = getf(opts, sub, field);
if ~present, return; end
checkTextValue([sub '.' field], v);
end

function checkTextValue(name, v)
if ~((ischar(v)) || (isstring(v) && isscalar(v)))
    error('ISA:ISAvalidateOpts:notText', 'opts.%s must be a char or scalar string; got %s.', ...
        name, describe(v));
end
end

function checkCellOfText(opts, sub, field)
[v, present] = getf(opts, sub, field);
if ~present, return; end
checkCellOfTextValue([sub '.' field], v);
end

function checkCellOfTextValue(name, v)
if ~(iscell(v) && all(cellfun(@(x) ischar(x) || (isstring(x) && isscalar(x)), v)))
    error('ISA:ISAvalidateOpts:notCellOfText', ...
        'opts.%s must be a cell array of char/string names; got %s.', name, describe(v));
end
end

function checkPositive(opts, sub, field, zeroAllowed)
if nargin < 4, zeroAllowed = false; end
[v, present] = getf(opts, sub, field);
if ~present, return; end
checkPositiveValue([sub '.' field], v, zeroAllowed);
end

function checkPositiveValue(name, v, zeroAllowed)
if nargin < 3, zeroAllowed = false; end
if ~(isnumeric(v) && isscalar(v) && isreal(v) && isfinite(v) && (v > 0 || (zeroAllowed && v == 0)))
    if zeroAllowed
        error('ISA:ISAvalidateOpts:notPositive', 'opts.%s must be a non-negative numeric scalar; got %s.', ...
            name, describe(v));
    else
        error('ISA:ISAvalidateOpts:notPositive', 'opts.%s must be a positive numeric scalar; got %s.', ...
            name, describe(v));
    end
end
end

function checkPosInt(opts, sub, field, zeroAllowed)
if nargin < 4, zeroAllowed = false; end
[v, present] = getf(opts, sub, field);
if ~present, return; end
checkPosIntValue([sub '.' field], v, zeroAllowed);
end

function checkPosIntValue(name, v, zeroAllowed)
if nargin < 3, zeroAllowed = false; end
checkPositiveValue(name, v, zeroAllowed);
if v ~= floor(v)
    error('ISA:ISAvalidateOpts:notInteger', 'opts.%s must be an integer; got %s.', name, describe(v));
end
end

function checkUnitRange(opts, sub, field)
[v, present] = getf(opts, sub, field);
if ~present, return; end
checkUnitRangeValue([sub '.' field], v);
end

function checkUnitRangeValue(name, v)
if ~(isnumeric(v) && isscalar(v) && isreal(v) && ~isnan(v) && v >= 0 && v <= 1)
    error('ISA:ISAvalidateOpts:notInUnitRange', ...
        'opts.%s must be a numeric scalar in [0, 1]; got %s.', name, describe(v));
end
end

function checkMember(opts, sub, field, validSet)
[v, present] = getf(opts, sub, field);
if ~present, return; end
checkMemberValue([sub '.' field], v, validSet);
end

function checkMemberValue(name, v, validSet)
if iscell(validSet) && ischar(validSet{1})
    ok = (ischar(v) || (isstring(v) && isscalar(v))) && any(strcmpi(char(v), validSet));
    optionsStr = ['{''' strjoin(validSet, ''', ''') '''}'];
else
    ok = isnumeric(v) && isscalar(v) && any(cellfun(@(x) isequal(v, x), validSet));
    optionsStr = ['{' strjoin(cellfun(@num2str, validSet, 'UniformOutput', false), ', ') '}'];
end
if ~ok
    error('ISA:ISAvalidateOpts:notMember', 'opts.%s must be one of %s; got %s.', ...
        name, optionsStr, describe(v));
end
end

function checkViewGroups(opts, sub, field)
[v, present] = getf(opts, sub, field);
if ~present, return; end
if ~iscell(v)
    error('ISA:ISAvalidateOpts:badViewGroups', ...
        'opts.pilot.viewGroups must be a cell array of algorithm index vectors; got %s.', describe(v));
end
for i = 1:numel(v)
    g = v{i};
    if ~(isnumeric(g) && isvector(g) && isreal(g) && all(isfinite(g)) && all(g > 0) && all(g == floor(g)))
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
