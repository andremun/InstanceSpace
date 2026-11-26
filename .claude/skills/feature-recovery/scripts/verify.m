function varargout = verify(featuresCsv, defsFile, varargin)
% verify  Check candidate feature definitions against the column values.
%
%   You write the definitions file; this checks whether each definition
%   actually reproduces its column, and when it does not, says what kind of
%   wrong it is.
%
%   verify('features.csv', 'definitions.m')
%   verify('features.csv', 'definitions.m', 'Columns', 'feature_n', ...
%          'Sample', 20)
%   r = verify(...)              % returns the verdicts instead of printing them
%
%   The definitions file is a function file that takes no arguments and
%   returns a struct with two fields:
%
%     load      handle, @(path) -> containers.Map from instance id to the
%               parsed instance, whatever type suits the problem
%     features  N-by-2 cell array, {column name, @(instance) -> scalar; ...}
%
%   load may take a second argument; when Sample is in use it is given the
%   ids actually being verified, as a cell array, so only those instances
%   need parsing. The file's own folder is on the path while verify runs, so
%   it can keep helpers in sibling files.
%
%   While searching for a reading, narrow the run rather than widening it:
%   Columns to the one column under study, Sample to 20 or 30 rows. A
%   subsample refutes a candidate; only the full run confirms one.
%
%   Options (name-value pairs):
%     Instances   path passed to load (default: the features CSV's folder)
%     IdCol       instance id column (default: the first column)
%     Quiet       only show what is not EXACT (default: false)
%     Columns     columns to verify, as a cell array or a comma-separated
%                 string (default: all of them)
%     Sample      verify on this many rows, drawn from the column extremes
%                 and spread across the source column (default: every row)
%     Seed        seed for the non-extreme part of the sample (default: 0)
%     Print       print the report (default: true when no output is requested)
%
%   Runs in MATLAB and in GNU Octave, with no toolbox.
%
%   Companion function for the feature-recovery skill.
%
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

% The definitions file's folder stays on the path for the run, so the
% definitions can call helpers kept beside them, and comes off again
% afterwards, however the run ends.
oldPath = path;
try
    [varargout{1:nargout}] = verifyRun(nargout, featuresCsv, defsFile, varargin{:});
catch err
    path(oldPath);
    rethrow(err);
end
path(oldPath);
end

function varargout = verifyRun(nout, featuresCsv, defsFile, varargin)
SHOW_UNDEFINED = 6;  % columns listed by name before the rest are summarised

p = inputParser;
p.FunctionName = 'verify';
addRequired(p, 'featuresCsv', @(x) ischar(x) || isstring(x));
addRequired(p, 'defsFile', @(x) ischar(x) || isstring(x));
addParameter(p, 'Instances', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'IdCol', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'Quiet', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'Columns', {}, @(x) ischar(x) || isstring(x) || iscell(x));
addParameter(p, 'Sample', [], @(x) isempty(x) || (isscalar(x) && isnumeric(x)));
addParameter(p, 'Seed', 0, @(x) isscalar(x) && isnumeric(x));
addParameter(p, 'Print', nout == 0, @(x) islogical(x) || isnumeric(x));
parse(p, featuresCsv, defsFile, varargin{:});
a = p.Results;
featuresCsv = char(a.featuresCsv);
say = makePrinter(a.Print);

% ---- the definitions --------------------------------------------------------
defs = loadDefinitions(char(a.defsFile));
defNames = defs.features(:, 1)';

% ---- the CSV ------------------------------------------------------------------
[header, cells] = readCells(featuresCsv);
if isempty(cells)
    error('verify:noRows', '%s has no data rows', featuresCsv);
end
idcol = char(a.IdCol);
if isempty(idcol), idcol = header{1}; end
idPos = find(strcmp(header, idcol), 1);
if isempty(idPos)
    error('verify:noIdCol', 'no column ''%s''; columns are %s', idcol, ...
        strjoin(header(1:min(8, end)), ', '));
end
isSrc = ismember(lower(header), {'source', 'sources'});
colPos = find(~isSrc & (1:numel(header)) ~= idPos);
cols = header(colPos);
want = a.Columns;
if isstring(want) && ~isscalar(want)
    want = cellstr(want);
elseif ischar(want) || isstring(want)
    want = strsplit(char(want), ',');
end
want = strtrim(cellstr(want));
want = want(:)';
want = want(~cellfun(@isempty, want));
if ~isempty(want)
    missing = want(~ismember(want, cols));
    if ~isempty(missing)
        error('verify:noColumn', 'no such column(s): %s; columns are %s', ...
            strjoin(missing, ', '), strjoin(cols(1:min(8, end)), ', '));
    end
    [~, k] = ismember(want, cols);
    colPos = colPos(k);
    cols = want;
end

% Whether each column varies over the *whole* CSV, recorded before any
% subsetting: a column that is flat only because the sample is small must
% not be reported as CONSTANT, which in this skill means unidentifiable.
varies = false(1, numel(cols));
for c = 1:numel(cols)
    varies(c) = numel(unique(cells(:, colPos(c)))) > 1;
end

% Sample against the columns actually being verified. The source column is
% excluded from scoring but is the best stratum there is: a missing branch
% is often one generator handled differently.
gPos = find(isSrc, 1);
fullCells = cells;
picked = [];
if ~isempty(a.Sample) && a.Sample < size(cells, 1)
    groups = {};
    if ~isempty(gPos), groups = cells(:, gPos); end
    [idx, spread] = chooseSample(cells, colPos, a.Sample, a.Seed, groups);
    cells = cells(idx, :);
    picked = [numel(idx), spread];
end

nsrc = 0;
if ~isempty(gPos), nsrc = numel(unique(fullCells(:, gPos))); end
ids = cells(:, idPos);
instPath = char(a.Instances);
if isempty(instPath), instPath = fileparts(absPath(featuresCsv)); end
narrowed = false;
if ~isempty(picked) && (nargin(defs.load) >= 2 || nargin(defs.load) < 0)
    objs = defs.load(instPath, ids);
    narrowed = true;
else
    objs = defs.load(instPath);
end
if ~isa(objs, 'containers.Map')
    error('verify:badLoad', 'load must return a containers.Map from id to instance');
end

keep = find(cellfun(@(x) isKey(objs, x), ids))';
if ~isempty(picked)
    covered = '';
    if ~isempty(gPos)
        covered = sprintf(', covering %d of %d sources', ...
            numel(unique(cells(:, gPos))), nsrc);
    end
    tail = '';
    if ~narrowed
        tail = ['; load parsed every instance -- give it an ids argument ' ...
                'to parse only these'];
    end
    say('sample    : %d rows, %d taken from the ends and quantiles of %d column(s)%s%s\n', ...
        picked(1), picked(2), numel(cols), covered, tail);
    say('            a sample refutes a reading; only the full run confirms one\n');
else
    say('alignment : %d/%d rows have instances; %d instances unused\n', ...
        numel(keep), numel(ids), objs.Count - numel(keep));
end
if isempty(keep)
    error('verify:noAlignment', 'nothing aligned -- check the id convention first');
end
if numel(keep) < numel(ids)
    say(['            verification denominator is %d, not %d; unmatched rows ' ...
         'can refute a definition but never confirm one\n'], numel(keep), numel(ids));
end
kids = ids(keep);
kobjs = values(objs, kids);

undefined = cols(~ismember(cols, defNames));
say('defined   : %d/%d columns\n\n', numel(cols) - numel(undefined), numel(cols));

% ---- score each column ----------------------------------------------------------
tally = struct('name', {}, 'count', {});
secsBy = zeros(1, 0);
secsNames = {};
res = struct('column', {}, 'verdict', {}, 'note', {}, 'failing', {});

for c = 1:numel(cols)
    col = cols{c};
    di = find(strcmp(defNames, col), 1);
    if isempty(di), continue; end
    s = strtrim(cells(keep, colPos(c)));
    [y, okNum] = parseNumbers(s);
    if ~all(okNum)
        bad = s(~okNum);
        fail(col, sprintf('%d non-numeric cell(s), e.g. %s', numel(bad), ...
            cellList(bad(1:min(3, end)))));
        continue;
    end
    if ~all(isfinite(y))
        % NaN and Inf parse cleanly and would then silently poison every
        % comparison, so stop here rather than report a meaningless verdict.
        fail(col, sprintf('%d non-finite cell(s) in the column', sum(~isfinite(y))));
        continue;
    end
    dec = cellfun(@decimals, s);
    isInt = all(y == round(y)) && all(dec == 0);
    tol = tolerances(s, dec, isInt);
    f = defs.features{di, 2};
    t0 = tic;
    try
        v = zeros(numel(kids), 1);
        for i = 1:numel(kids)
            out = f(kobjs{i});
            if ~isnumeric(out) && ~islogical(out) || numel(out) ~= 1
                error('verify:notScalar', 'definition returned a %s of size %s', ...
                    class(out), mat2str(size(out)));
            end
            v(i) = double(out);
        end
    catch err
        fail(col, sprintf('%s: %s', errName(err), err.message));
        continue;
    end
    secsBy(end + 1) = toc(t0); %#ok<AGROW>
    secsNames{end + 1} = col; %#ok<AGROW>
    if isscalar(unique(round12(y)))
        if ~isempty(picked) && varies(c)
            say('%-52s FLAT        constant in this sample but not in the column; raise Sample to score it\n', col);
            tally = bump(tally, 'FLAT');
            res(end + 1) = struct('column', col, 'verdict', 'FLAT', 'note', '', 'failing', 0); %#ok<AGROW>
            continue;
        end
        say('%-52s CONSTANT    column never varies; any constant matches it\n', col);
        tally = bump(tally, 'CONSTANT');
        res(end + 1) = struct('column', col, 'verdict', 'CONSTANT', 'note', '', 'failing', 0); %#ok<AGROW>
        continue;
    end
    [verdict, note, bad] = diagnose(y, v, tol);
    tally = bump(tally, verdict);
    res(end + 1) = struct('column', col, 'verdict', verdict, 'note', note, ...
        'failing', sum(bad)); %#ok<AGROW>
    if ismember(verdict, {'EXACT', 'EXACT~', 'EXACT~f32'}) && a.Quiet
        continue;
    end
    secs = '';
    if secsBy(end) >= 1.0, secs = sprintf('  [%.1fs]', secsBy(end)); end
    say('%-52s %-11s %s%s\n', col, verdict, note, secs);
    if ~isempty(bad) && ~strcmp(verdict, 'NO MATCH')
        line = structured(bad, kids, y, v);
        if ~isempty(line), say('%s\n', line); end
    end
end

[~, order] = sort({tally.name});
parts = arrayfun(@(t) sprintf('%s: %d', t.name, t.count), tally(order), ...
    'UniformOutput', false);
say('\n%s\n', strjoin(parts, '  '));
if ~isempty(undefined)
    more = '';
    if numel(undefined) > SHOW_UNDEFINED, more = ' ...'; end
    say('no definition supplied for %d: %s%s\n', numel(undefined), ...
        cellList(undefined(1:min(6, end))), more);
end

% Which column to go and make cheaper -- a guess about where the time went
% is usually wrong, and the answer decides whether Sample is worth it.
total = sum(secsBy);
if total >= 1.0
    [t, o] = sort(secsBy, 'descend');
    o = o(t >= 0.05 * total);
    o = o(1:min(3, end));
    slow = arrayfun(@(k) sprintf('%s %.1fs', secsNames{k}, secsBy(k)), o, ...
        'UniformOutput', false);
    say('time      : %.1fs over %d instances (%.1f ms/instance); slowest %s\n', ...
        total, numel(keep), total / max(1, numel(keep)) * 1e3, strjoin(slow, ', '));
end

if nout > 0
    r.columns = res;
    r.tally = tally(order);
    r.aligned = numel(keep);
    r.rows = numel(ids);
    r.undefined = undefined;
    varargout{1} = r;
end

    function fail(col, msg)
        % One column's problem must not end the run: report it and carry on.
        say('%-52s ERROR       %s\n', col, msg);
        tally = bump(tally, 'ERROR');
        res(end + 1) = struct('column', col, 'verdict', 'ERROR', 'note', msg, 'failing', 0);
    end
end

% =================================================================== tolerance ==

function d = decimals(s)
% Decimal places a CSV cell was printed with, or NaN for scientific notation.
if any(s == 'e' | s == 'E')
    d = NaN;
else
    k = find(s == '.', 1);
    if isempty(k), d = 0; else, d = numel(s) - k; end
end
end

function tol = tolerances(s, dec, isInt)
% The rounding tolerance for each cell of a column.
%   Per value, not per column: a CSV at fixed significant figures carries
%   less absolute precision at larger magnitudes, so a column-wide bound
%   wrongly rejects the big rows. Slack keeps an exact half-ulp tie passing.
if isInt
    tol = repmat(1e-9, numel(s), 1);
    return;
end
x = abs(str2double(s(:)));
rounded = 0.5 * 10 .^ -dec(:) * (1 + 1e-6);
sci = isnan(dec(:));
rounded(sci) = 1e-9 * max(1, x(sci));
% A CSV written at full double precision prints the exact double, so the
% rounding bound collapses below one ulp. Never demand better than float
% reproduction: recomputing reorders operations and moves the last bits.
tol = max(rounded, 8 * eps * max(1, x));
end

function t = loose(y)
% Tolerance at single precision, for telling rounding from a wrong reading.
%   A column stored as float32 before it was printed carries about seven
%   significant digits however many decimals the CSV shows, so a definition
%   that is right still lands this far out. Storage width cannot be read
%   back once a CSV has rounded to fixed decimal places, so every column is
%   scored at both bounds and told apart by which one it meets.
t = 1e-6 * max(1, abs(y));
end

% ================================================================= diagnostics ==

function note = partialNote(hit, okf, n)
% Describe a partial match at both tolerances. The failures are diagnosed at
% the looser one, since rows that differ only by rounding are not the ones
% to go and characterise.
note = sprintf('matches %d/%d', hit, n);
if sum(okf) > hit
    note = sprintf(['%s; %d/%d to ~7 significant digits, so %d are real ' ...
        'disagreements and the rest is rounding'], note, sum(okf), n, n - sum(okf));
end
end

function [verdict, note, bad] = diagnose(y, v, tol)
% Classify how a definition's values compare with the column. Each verdict
% points somewhere specific, so none of them is 'try again'.
NOISE_REL = 1e-10;          % last-bits disagreement between two exact computations
EXACT_FIT = 1e-12;          % correlation, slope and offset indistinguishable from exact
UNIT_SLOPE = 1e-9;          % a fitted slope this close to 1 is no rescaling at all
NEAR_PROPORTIONAL = 0.999999;
SAME_ORDER = 0.9999;

ok = abs(v - y) <= tol;
n = numel(y);
hit = sum(ok);
note = '';
bad = [];
if hit == n
    verdict = 'EXACT';
    return;
end
if ~all(isfinite(v))
    verdict = 'ERROR';
    note = 'definition returned non-finite values';
    return;
end
if isscalar(unique(round12(v))) && ~isscalar(unique(round12(y)))
    verdict = 'NO MATCH';
    note = 'definition is constant; the column is not';
    bad = ~ok;
    return;
end

% Same quantity, last bits apart: iterative solvers and long accumulations
% cannot reproduce bit-for-bit, so this is a pass with a note rather than a
% miss -- kept distinct from a genuine scale error.
scale = max(1, abs(y));
if max(abs(v - y) ./ scale) < NOISE_REL
    verdict = 'EXACT~';
    note = sprintf('within numerical noise (max |Δ| = %s)', fmtG(max(abs(v - y)), 2));
    return;
end

% Every row agrees to single precision. A column printed with more decimals
% than it was stored with cannot do better, and neither can a long
% accumulation, so the reading is right and the residue is rounding.
okf = abs(v - y) <= loose(y);
if all(okf)
    verdict = 'EXACT~f32';
    note = sprintf(['agrees to ~7 significant digits (max rel = %s); the column ' ...
        'is not stored more precisely than this'], fmtE(max(abs(v - y) ./ scale), 1));
    return;
end

% A reading exact on most rows is missing a branch, not a transform, however
% well it correlates -- a few far-out rows drag the fit below. Count only
% the non-zero rows: a wrong scale matches every zero for free.
nz = abs(y) > tol;
if any(nz) && sum(okf & nz) * 2 > sum(nz)
    verdict = 'PARTIAL';
    note = partialNote(hit, okf, n);
    bad = ~okf;
    return;
end

% affine: right quantity, wrong scale -- the fitted slope names the gap
if std(v, 1) > 0
    ab = [v, ones(n, 1)] \ y;
    sa = ab(1);
    sb = ab(2);
    cc = corrcoef(v, y);
    r = cc(1, 2);
    if r > 1 - EXACT_FIT
        hint = sprintf('scale by %s', fmtG(sa, 10));
        if abs(sb) > EXACT_FIT
            hint = sprintf('%s, offset %s', hint, fmtG(sb, 6, true));
        end
        if abs(sa - 1) > UNIT_SLOPE || abs(sb) > EXACT_FIT
            verdict = 'AFFINE';
            note = sprintf('exactly proportional -- %s', hint);
            bad = ~ok;
            return;
        end
    end
    if r > NEAR_PROPORTIONAL
        verdict = 'AFFINE';
        note = sprintf('near-proportional (r=%.9f), scale %s', r, fmtG(sa, 8));
        bad = ~ok;
        return;
    end
    rho = spearman(v, y);
    if rho > SAME_ORDER
        verdict = 'MONOTONE';
        note = sprintf('same ordering (rho=%.6f); unknown transform', rho);
        bad = ~ok;
        return;
    end
end
if hit
    verdict = 'PARTIAL';
    note = partialNote(hit, okf, n);
    bad = ~okf;
elseif any(okf)
    verdict = 'PARTIAL';
    note = sprintf('matches %d/%d to ~7 significant digits', sum(okf), n);
    bad = ~okf;
else
    verdict = 'NO MATCH';
    note = 'no relationship';
    bad = ~ok;
end
end

function line = structured(bad, ids, y, v)
% Summarise the failing rows of a column. Scattered failures are edge cases;
% a predictable set is a missing branch. Reported with the failing ids so
% the pattern can be characterised -- and weigh any claimed pattern against
% how common it is overall.
k = sum(bad);
line = '';
if k == 0, return; end
w = find(bad);
ex = ids(w(1:min(4, end)));
d = abs(v(bad) - y(bad));
line = sprintf('      %d failing, e.g. %s  |Δ| median %s max %s', k, cellList(ex), ...
    fmtG(median(d), 3), fmtG(max(d), 3));
end

% ====================================================================== sample ==

function out = bisectOrder(n)
% Positions 1..n ordered ends-first, then repeatedly bisecting. The first k
% of it is a sensible spread for any k, so raising Sample extends the
% sample instead of replacing it, and a row that refuted a candidate at 20
% is still there at 50.
if n <= 0, out = zeros(1, 0); return; end
if n == 1, out = 1; return; end
out = [1, n];
seen = false(1, n);
seen([1, n]) = true;
q = [1, n];
head = 1;
while head <= size(q, 1)
    lo = q(head, 1);
    hi = q(head, 2);
    head = head + 1;
    if hi - lo < 2, continue; end
    mid = floor((lo - 1 + hi - 1) / 2) + 1;
    if ~seen(mid)
        out(end + 1) = mid; %#ok<AGROW>
        seen(mid) = true;
    end
    q(end + 1, :) = [lo, mid]; %#ok<AGROW>
    q(end + 1, :) = [mid, hi]; %#ok<AGROW>
end
end

function lane = valueLane(cells, idx, pos)
% Rows of idx ordered ends-first by the column at pos, or [] if it is not numeric.
[vals, ok] = parseNumbers(cells(idx, pos));
if ~all(ok), lane = []; return; end
[~, o] = sort(vals);
order = idx(o);
lane = order(bisectOrder(numel(order)));
lane = lane(:)';
end

function [picks, spread] = chooseSample(cells, colPos, n, seed, groups)
% Choose the rows a Sample run verifies: n rows spread across the value range
% of each column under study, and across the sources if the CSV names one
% per row.
%
% Not a uniform draw. What a subsample has to catch is a missing branch, and
% those sit at the ends of a column's range -- the degenerate instance, the
% one row exactly on a boundary, the largest instance -- which random
% sampling reaches with low probability precisely because they are rare.
%
% Two strata, in priority order. The global extremes of each column go in
% first, since a boundary convention shows up there and nowhere else. Then
% the picks are spread across the sources rather than left to cluster in
% one, because the other common missing branch is a whole generator handled
% differently. Within a source, rows are again taken ends-first.
%
% Deterministic, so two candidate readings are compared on the same rows.
m = size(cells, 1);
if n >= m
    picks = 1:m;
    spread = 0;
    return;
end
picks = zeros(1, 0);
seen = false(1, m);

colLanes = {};
for c = 1:numel(colPos)
    lane = valueLane(cells, 1:m, colPos(c));
    if ~isempty(lane), colLanes{end + 1} = lane; end %#ok<AGROW>
end

% The extremes first -- but capped, so that on a wide CSV they cannot eat
% the whole budget and crowd out the sources.
for c = 1:numel(colLanes)
    if numel(picks) >= max(2, floor(n / 2)), break; end
    lane = colLanes{c};
    for i = lane(1:min(2, end))
        take(i);
    end
end

% One row per source before any source gets a second.
srcLanes = {};
if ~isempty(groups)
    [g, ~, gi] = unique(groups);
    for k = 1:numel(g)
        members = find(gi == k)';
        lane = [];
        if ~isempty(colPos), lane = valueLane(cells, members, colPos(1)); end
        if isempty(lane), lane = members; end
        srcLanes{end + 1} = lane; %#ok<AGROW>
    end
end

laneSets = {srcLanes, colLanes};
for s = 1:2
    lanes = laneSets{s};
    depth = max([0, cellfun(@numel, lanes)]);
    for tier = 1:depth
        if numel(picks) >= n, break; end
        for k = 1:numel(lanes)
            if tier <= numel(lanes{k}), take(lanes{k}(tier)); end
        end
    end
end
spread = numel(picks);

rest = find(~seen);
need = max(0, n - numel(picks));
if need > 0
    prev = rng;
    restore = onCleanup(@() rng(prev));
    rng(seed, 'twister');
    rest = rest(randperm(numel(rest)));
    clear restore
    picks = [picks, rest(1:need)];
end
picks = sort(picks);

    function take(i)
        if ~seen(i) && numel(picks) < n
            picks(end + 1) = i;
            seen(i) = true;
        end
    end
end

% ===================================================================== helpers ==

function defs = loadDefinitions(file)
% Run the definitions file and check what it returned.
file = absPath(file);
[folder, name, ext] = fileparts(file);
if ~strcmp(ext, '.m') || exist(file, 'file') ~= 2
    error('verify:noDefs', 'cannot find the definitions file %s', file);
end
% Its own folder goes on the path first, so the definitions can split
% helpers -- an expensive kernel, a wrapper round a compiled one -- into
% sibling files instead of being one file by necessity. It is called from
% inside that folder, because a file of the same name in the current folder
% would otherwise take precedence over it.
addpath(folder);
here = cd(folder);
back = onCleanup(@() cd(here));
defs = feval(name);
clear back
if ~isstruct(defs) || ~isfield(defs, 'load') || ~isfield(defs, 'features')
    error('verify:badDefs', '%s must return a struct with fields load and features', file);
end
if ~isa(defs.load, 'function_handle')
    error('verify:badDefs', '%s: load must be a function handle', file);
end
f = defs.features;
if ~iscell(f) || size(f, 2) ~= 2 || ~all(cellfun(@ischar, f(:, 1))) ...
        || ~all(cellfun(@(h) isa(h, 'function_handle'), f(:, 2)))
    error('verify:badDefs', ['%s: features must be an N-by-2 cell array of ' ...
        '{column name, function handle}'], file);
end
end

function [header, cells] = readCells(path)
% Every cell of a CSV as text, since the printed digits set the tolerance.
% The delimiter is whichever of comma, semicolon or tab the header holds most of.
if exist(path, 'file') ~= 2
    error('verify:cannotOpen', 'cannot open %s', path);
end
txt = fileread(path);
if ~isempty(txt) && double(txt(1)) == 65279, txt = txt(2:end); end  % BOM
lines = regexp(txt, '\r?\n', 'split');
lines = lines(~cellfun(@(l) isempty(strtrim(l)), lines));
if isempty(lines)
    error('verify:noRows', '%s is empty', path);
end
delims = {',', ';', sprintf('\t')};
counts = cellfun(@(d) sum(lines{1} == d), delims);
[top, best] = max(counts);
delim = ',';
if top > 0, delim = delims{best}; end
header = splitLine(lines{1}, delim);
cells = cell(numel(lines) - 1, numel(header));
for i = 2:numel(lines)
    f = splitLine(lines{i}, delim);
    f(end + 1:numel(header)) = {''};
    cells(i - 1, :) = f(1:numel(header));
end
end

function f = splitLine(line, delim)
% Split one CSV line, honouring double-quoted fields.
if ~any(line == '"')
    f = strsplit(line, delim, 'CollapseDelimiters', false);
    return;
end
f = {};
buf = '';
inq = false;
i = 1;
while i <= numel(line)
    ch = line(i);
    if inq
        if ch == '"' && i < numel(line) && line(i + 1) == '"'
            buf(end + 1) = '"'; %#ok<AGROW>
            i = i + 1;
        elseif ch == '"'
            inq = false;
        else
            buf(end + 1) = ch; %#ok<AGROW>
        end
    elseif ch == '"'
        inq = true;
    elseif ch == delim
        f{end + 1} = buf; %#ok<AGROW>
        buf = '';
    else
        buf(end + 1) = ch; %#ok<AGROW>
    end
    i = i + 1;
end
f{end + 1} = buf;
end

function [x, ok] = parseNumbers(s)
% Parse text cells as numbers. NaN and Inf spelled out count as numbers,
% which is what lets a non-finite column be reported as such.
x = str2double(s(:));
ok = ~isnan(x) & imag(x) == 0;
x = real(x);
lit = regexpi(strtrim(s(:)), '^[+-]?(nan|inf|infinity)$', 'once');
ok = ok | ~cellfun(@isempty, lit);
end

function y = round12(x)
y = round(x * 1e12) / 1e12;
end

function rho = spearman(a, b)
% Pearson correlation of average ranks.
ra = avgRank(a);
rb = avgRank(b);
cc = corrcoef(ra, rb);
rho = cc(1, 2);
end

function r = avgRank(x)
[s, o] = sort(x(:));
n = numel(s);
r = zeros(n, 1);
i = 1;
while i <= n
    j = i;
    while j < n && s(j + 1) == s(i), j = j + 1; end
    r(o(i:j)) = (i + j) / 2;
    i = j + 1;
end
end

function tally = bump(tally, name)
k = find(strcmp({tally.name}, name), 1);
if isempty(k)
    tally(end + 1) = struct('name', name, 'count', 1);
else
    tally(k).count = tally(k).count + 1;
end
end

function s = cellList(c)
% A cell array of text written as a list: {'a', 'b'}.
s = ['{' strjoin(cellfun(@(x) ['''' x ''''], c(:)', 'UniformOutput', false), ', ') '}'];
end

function s = errName(err)
s = err.identifier;
if isempty(s), s = 'error'; end
end

function s = fmtG(x, sig, signed)
% %g at sig significant digits, with a two-digit exponent.
if nargin < 3, signed = false; end
if signed, f = sprintf('%%+.%dg', sig); else, f = sprintf('%%.%dg', sig); end
s = regexprep(sprintf(f, x), 'e([+-])0*(\d\d)', 'e$1$2');
end

function s = fmtE(x, digits)
s = regexprep(sprintf(sprintf('%%.%de', digits), x), 'e([+-])0*(\d\d)', 'e$1$2');
end

function p = absPath(p)
% The absolute form of a path, for files that exist.
if ~(numel(p) >= 1 && (p(1) == '/' || p(1) == '\' || (numel(p) >= 2 && p(2) == ':')))
    p = fullfile(pwd, p);
end
end

function say = makePrinter(enabled)
if enabled
    say = @(varargin) fprintf(varargin{:});
else
    say = @(varargin) [];
end
end
