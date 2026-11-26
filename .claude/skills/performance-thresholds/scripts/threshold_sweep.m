function varargout = threshold_sweep(csvPath, varargin)
% threshold_sweep  Performance-only diagnostics for ISA good-performance thresholds.
%
%   Answers one question: if you set opts.perf to these values, what labels
%   will PRELIM actually produce, and are they any good?
%
%   r = threshold_sweep('metadata.csv', 'MaxPerf', false, 'AbsPerf', true, ...
%                       'Epsilon', [1 10 60 600 3599]);
%   threshold_sweep('metadata.csv', 'MaxPerf', false, 'AbsPerf', true, ...
%                   'Epsilon', [1 10 60 600 3599])   % no output: prints the report
%
%   Reads a standard ISA metadata CSV and takes its algo_* columns, plus a
%   Source column when one is present.
%
%   Options (name-value pairs):
%     MaxPerf     true = higher is better (default false)
%     AbsPerf     true = fixed threshold in the metric's units, false =
%                 relative to the per-instance best (default false)
%     Epsilon     candidate thresholds (default: a fallback sweep grid)
%     PlateauTol  good-rate tolerance behind plateauWidth (default 0.05)
%     GridPoints  points on the dense good-rate curve (default 200)
%     Print       print the report (default: true when no output is requested)
%
%   Output r, a struct:
%     data          threshold-independent degeneracy report -- read first
%     perAlgo       table of per-algorithm min, max, median, nanFrac, range
%     sweep         table, one row per candidate threshold
%     curve         table, the dense good-rate curve
%     opts          the options used
%     labels        cell array, the Ybin matrix for each candidate
%     uniqueByAlgo  cell array, each algorithm's unique-good share per candidate
%
%   The labelling is a deliberate copy of core/PRELIM.m lines 89-120 rather
%   than an improvement on it: the job is to predict what PRELIM will
%   produce, quirks included, so a cleaner reimplementation would be a worse
%   tool. Verify with:
%
%     [~, ~, out] = PRELIM(X, Y, opts);
%     isequal(r.labels{k}, out.Ybin)
%
%   Runs in MATLAB and in GNU Octave. Octave needs the datatypes package
%   loaded for the result tables, as the toolkit itself does.
%
%   Companion function for the performance-thresholds skill.
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

% Good-rates outside [RATE_LOW, RATE_HIGH] count as saturated or near-degenerate.
RATE_LOW = 0.05;
RATE_HIGH = 0.95;
% PRELIM's cut for too many zero best-scores, and the share of values a point
% mass needs before it is treated as a sentinel.
ZERO_BEST_FRAC = 0.05;
POINT_MASS_FRAC = 0.05;
% An algorithm uniquely good on less than this share of instances is a passenger.
PASSENGER_FRAC = 0.01;
% Fewest algorithms for a screen score.
MIN_SCREEN_ALGOS = 3;

p = inputParser;
p.FunctionName = 'threshold_sweep';
addRequired(p, 'csvPath', @(x) ischar(x) || isstring(x));
addParameter(p, 'MaxPerf', false, @(x) isscalar(x) && (islogical(x) || isnumeric(x)));
addParameter(p, 'AbsPerf', false, @(x) isscalar(x) && (islogical(x) || isnumeric(x)));
addParameter(p, 'Epsilon', [], @isnumeric);
addParameter(p, 'PlateauTol', 0.05, @(x) isscalar(x) && isnumeric(x));
addParameter(p, 'GridPoints', 200, @(x) isscalar(x) && isnumeric(x));
addParameter(p, 'Print', nargout == 0, @(x) isscalar(x) && (islogical(x) || isnumeric(x)));
parse(p, csvPath, varargin{:});
o = p.Results;
maxPerf = logical(o.MaxPerf);
absPerf = logical(o.AbsPerf);

[Y, algoLabels, groups] = readAlgoColumns(char(o.csvPath));
[ninst, nalgos] = size(Y);
finite = Y(isfinite(Y));
ybest = bestPerInstance(Y, maxPerf);

uniqVals = unique(finite);
if numel(uniqVals) > 1
    resolution = min(diff(uniqVals));
else
    resolution = NaN;
end
if ~isempty(finite)
    [vals, ~, ic] = unique(finite);
    counts = accumarray(ic(:), 1);
    [cmax, imax] = max(counts);
    pmVal = vals(imax);
    pmFrac = cmax / numel(finite);
else
    pmVal = NaN;
    pmFrac = 0;
end

% A sentinel (timeout, ceiling) needs a rung strictly inside it, and how far
% inside is not free: any genuine value between the rung and the sentinel is
% relabelled as a censored run. These are the nearest genuine values on each
% side of the point mass, so the caller can see the interval a rung may
% safely occupy instead of guessing at cutoff - 1.
below = uniqVals(uniqVals < pmVal);
above = uniqVals(uniqVals > pmVal);
pmBelow = NaN;
pmAbove = NaN;
if ~isempty(below), pmBelow = max(below); end
if ~isempty(above), pmAbove = min(above); end

% Ties matter whenever epsilon = 0 is a live candidate: at that rung every
% tied-best algorithm is labelled good, so the tie rate is the entire
% deviation from a clean winner-take-all split. Mirrors the multiple-best
% statistic PRELIM reports.
tieFrac = mean(sum(bsxfun(@eq, Y, ybest), 2) > 1);

data = struct();
data.nInstances = ninst;
data.nAlgos = nalgos;
data.algoLabels = algoLabels;
data.globalMin = extremum(finite, @min);
data.globalMax = extremum(finite, @max);
data.nanFrac = mean(isnan(Y(:)));
data.ybestZeroFrac = mean(ybest == 0);
data.tieFrac = tieFrac;
data.hasNegative = any(finite < 0);
data.crossesZero = any(finite < 0) && any(finite > 0);
data.pointMassValue = pmVal;
data.pointMassFrac = pmFrac;
data.pointMassNeighbourBelow = pmBelow;
data.pointMassNeighbourAbove = pmAbove;
data.resolution = resolution;
data.nDistinct = numel(uniqVals);

% When the point mass sits at an extreme it is almost always a sentinel, and
% the rung has to clear it without swallowing genuine runs. Name the safe
% interval explicitly -- the failure it prevents is silent.
% A sentinel that must be labelled bad has to sit on the bad side of the
% rung, which pins it to the max under minimising and to the min under
% maximising; those are the only two arrangements a rung can resolve. A mass
% on the *good* side cannot be cleared by any epsilon, so say that instead of
% returning nothing.
if pmFrac >= POINT_MASS_FRAC && isfinite(pmVal)
    atMax = pmVal == data.globalMax;
    atMin = pmVal == data.globalMin;
    if atMax && ~maxPerf
        data.safeRungInterval = [pmBelow, pmVal];  % minimising: (gap, sentinel)
    elseif atMin && maxPerf
        data.safeRungInterval = [pmVal, pmAbove];  % maximising: (sentinel, gap)
    elseif atMax || atMin
        data.safeRungInterval = 'goodSide';
    else
        data.safeRungInterval = [];
    end
end

% Relative mode is meaningless on a scale that hits or crosses zero. PRELIM
% substitutes eps and warns at the same 5% cut; warn here too so the caller
% sees it before proposing a mode.
if ~absPerf && data.ybestZeroFrac > ZERO_BEST_FRAC
    warning('threshold_sweep:manyZeroBest', ['%.1f%% of instances have a ' ...
        'best-algorithm score of exactly zero. Relative mode is degenerate ' ...
        'on this data -- use AbsPerf = true.'], 100 * data.ybestZeroFrac);
end
if ~absPerf && data.crossesZero
    warning('threshold_sweep:crossesZero', ['Performance values cross zero, ' ...
        'so the relative ratio has no consistent meaning. Use AbsPerf = true.']);
end

amin = nan(nalgos, 1);
amax = nan(nalgos, 1);
amed = nan(nalgos, 1);
for k = 1:nalgos
    col = Y(:, k);
    col = col(~isnan(col));
    if ~isempty(col)
        amin(k) = min(col);
        amax(k) = max(col);
        amed(k) = median(col);
    end
end
perAlgo = table(algoLabels(:), amin, amax, amed, mean(isnan(Y), 1)', amax - amin, ...
    'VariableNames', {'algorithm', 'min', 'max', 'median', 'nanFrac', 'range'});

% Dense good-rate curve, computed independently of the reported rungs so
% plateau width and local slope are not limited by a coarse ladder.
if absPerf
    lo = data.globalMin;
    hi = data.globalMax;
    if ~isfinite(lo) || ~isfinite(hi) || lo == hi
        lo = 0;
        hi = 1;
    end
    pad = 0.02 * (hi - lo);
    grid = linspace(lo - pad, hi + pad, o.GridPoints);
else
    grid = linspace(0, 1, o.GridPoints);
end
gRate = zeros(size(grid));
for g = 1:numel(grid)
    lab = binarise(Y, maxPerf, absPerf, grid(g));
    gRate(g) = mean(lab(:));
end
curve = table(grid(:), gRate(:), 'VariableNames', {'epsilon', 'goodRate'});

% Candidate rungs. The default is an explicit fallback for sweeping, NOT a
% proposal -- the ladder is supposed to come from the domain (skill section 3).
if isempty(o.Epsilon)
    if absPerf
        rungs = ladder125(data.globalMin, data.globalMax);
    else
        rungs = [0 0.01 0.02 0.05 0.1 0.2 0.5];
    end
else
    rungs = o.Epsilon;
end
rungs = sort(double(rungs(:)'));

% Absolute rungs are in the metric's own units; relative rungs are fractions
% of the per-instance best. Passing a relative ladder (0.01, 0.05, ...) while
% asking for absolute mode is the common slip, and it fails silently -- every
% rung below the data's minimum labels everything bad, every rung above its
% maximum labels everything good.
if absPerf && isfinite(data.globalMin)
    stray = rungs(rungs < data.globalMin | rungs > data.globalMax);
    if ~isempty(stray)
        warning('threshold_sweep:strayRungs', ['%d of %d absolute rungs fall ' ...
            'outside the data range [%.6g, %.6g]: %s. Absolute rungs are in the ' ...
            'metric''s own units. If you meant fractions of the per-instance ' ...
            'best, use AbsPerf = false.'], numel(stray), numel(rungs), ...
            data.globalMin, data.globalMax, strtrim(sprintf('%.6g ', stray)));
    end
end

[gidx, ngroups, ordered] = groupIndex(groups);
nr = numel(rungs);
allLabels = cell(1, nr);
uniqByAlgo = cell(1, nr);
cols = struct('epsilon', rungs(:), 'goodRate', nan(nr, 1), 'algoRateMin', nan(nr, 1), ...
    'algoRateMax', nan(nr, 1), 'nDegenerateAlgos', nan(nr, 1), 'fracAllGood', nan(nr, 1), ...
    'fracAllBad', nan(nr, 1), 'meanJaccard', nan(nr, 1), 'uniqueBestFrac', nan(nr, 1), ...
    'groupSpread', nan(nr, 1), 'groupRho', nan(nr, 1), 'nNearDegenerate', nan(nr, 1), ...
    'nPassengers', nan(nr, 1), 'dGoodRate', nan(nr, 1), 'plateauWidth', nan(nr, 1), ...
    'screenScore', nan(nr, 1));

for ri = 1:nr
    e = rungs(ri);
    lab = binarise(Y, maxPerf, absPerf, e);
    allLabels{ri} = lab;
    goodRate = mean(lab(:));
    algoRate = mean(lab, 1);
    ngood = sum(lab, 2);
    fracAllGood = mean(ngood == nalgos);
    fracAllBad = mean(ngood == 0);
    jac = meanPairwiseJaccard(lab);

    % Screening components. Balance peaks at a 50% good-rate;
    % informativeness is the share of instances that discriminate between at
    % least two algorithms; complementarity rewards non-overlapping sets.
    balance = 1 - abs(2 * goodRate - 1);
    inform = 1 - fracAllGood - fracAllBad;
    if isnan(jac)
        compl = 0;
    else
        compl = 1 - jac;
    end
    screen = nthroot(max(balance, 0) * max(inform, 0) * max(compl, 0), 3);
    % With two algorithms fracAllBad is structurally zero and meanJaccard
    % equals fracAllGood, so the score collapses to a monotone function of
    % epsilon peaking at 0 -- it measures disagreement, not discrimination.
    if nalgos < MIN_SCREEN_ALGOS
        screen = NaN;
    end

    % Per-algorithm unique contribution: instances this algorithm is the only
    % one to get right. An algorithm with a decent good-rate but near-zero
    % unique contribution is a passenger.
    uniq = zeros(1, nalgos);
    for a = 1:nalgos
        others = lab(:, [1:a-1, a+1:nalgos]);
        uniq(a) = mean(lab(:, a) & ~any(others, 2));
    end
    uniqByAlgo{ri} = uniq;

    % Degeneracy guard: does this rung sort instance classes rather than
    % algorithms? Rejection criterion only -- never an objective.
    gSpread = NaN;
    gRho = NaN;
    if ngroups > 1
        gr = zeros(1, ngroups);
        for g = 1:ngroups
            block = lab(gidx == g, :);
            gr(g) = mean(block(:));
        end
        gSpread = max(gr) - min(gr);
        if ordered
            gRho = spearmanRho(1:ngroups, gr);
        end
    end

    [slope, width] = plateauStats(grid, gRate, e, goodRate, o.PlateauTol, RATE_LOW, RATE_HIGH);
    cols.goodRate(ri) = goodRate;
    cols.algoRateMin(ri) = min(algoRate);
    cols.algoRateMax(ri) = max(algoRate);
    cols.nDegenerateAlgos(ri) = sum(algoRate == 0 | algoRate == 1);
    cols.fracAllGood(ri) = fracAllGood;
    cols.fracAllBad(ri) = fracAllBad;
    cols.meanJaccard(ri) = jac;
    cols.uniqueBestFrac(ri) = mean(ngood == 1);
    cols.groupSpread(ri) = gSpread;
    cols.groupRho(ri) = gRho;
    % Exact 0/1 misses an algorithm good on 91% or 5% of instances, which is
    % still a badly unbalanced target for PYTHIA. Informational, not
    % disqualifying.
    cols.nNearDegenerate(ri) = sum(algoRate <= RATE_LOW | algoRate >= RATE_HIGH);
    cols.nPassengers(ri) = sum(uniq < PASSENGER_FRAC);
    cols.dGoodRate(ri) = slope;
    cols.plateauWidth(ri) = width;
    cols.screenScore(ri) = screen;
end

names = fieldnames(cols);
sweep = table();
for k = 1:numel(names)
    sweep.(names{k}) = cols.(names{k});
end

r = struct();
r.data = data;
r.perAlgo = perAlgo;
r.sweep = sweep;
r.curve = curve;
r.opts = struct('MaxPerf', maxPerf, 'AbsPerf', absPerf, ...
    'PlateauTol', o.PlateauTol, 'GridPoints', o.GridPoints);
r.labels = allLabels;
r.uniqueByAlgo = uniqByAlgo;

if o.Print
    printReport(r);
end
if nargout > 0
    varargout{1} = r;
end
end

% =========================================================================
% Core labelling -- a deliberate copy of core/PRELIM.m lines 89-120
% =========================================================================
function ybest = bestPerInstance(Y, maxPerf)
% bestPerInstance  Per-instance best score, matching PRELIM's NaN handling.
Yaux = Y;
if maxPerf
    Yaux(isnan(Yaux)) = -Inf;
    ybest = max(Yaux, [], 2);
else
    Yaux(isnan(Yaux)) = Inf;
    ybest = min(Yaux, [], 2);
end
end

function Ybin = binarise(Y, maxPerf, absPerf, epsilon)
% binarise  Good/bad labels, replicating core/PRELIM.m exactly.
%   The relative branch labels from Y with its zeros replaced by eps and its
%   NaNs left in place, divided by the eps-substituted best, as PRELIM does.
if maxPerf
    Yaux = Y;
    Yaux(isnan(Yaux)) = -Inf;
    Ybest = max(Yaux, [], 2);
    if absPerf
        Ybin = Yaux >= epsilon;
    else
        Ybest(Ybest == 0) = eps;
        Y(Y == 0) = eps;
        Ybin = (1 - bsxfun(@rdivide, Y, Ybest)) <= epsilon;
    end
else
    Yaux = Y;
    Yaux(isnan(Yaux)) = Inf;
    Ybest = min(Yaux, [], 2);
    if absPerf
        Ybin = Yaux <= epsilon;
    else
        Ybest(Ybest == 0) = eps;
        Y(Y == 0) = eps;
        Ybin = (bsxfun(@rdivide, Y, Ybest) - 1) <= epsilon;
    end
end
end

% =========================================================================
% Helpers
% =========================================================================
function [Y, labels, groups] = readAlgoColumns(csvPath)
% readAlgoColumns  Pull the algo_* columns (and any Source column) from an ISA CSV.
[names, cols] = readCsv(csvPath);
isalgo = strncmpi(names, 'algo_', 5);
if ~any(isalgo)
    error('threshold_sweep:noAlgoColumns', ...
        'No columns starting with "algo_" found in %s.', csvPath);
end
algoNames = names(isalgo);
if ~all(cellfun(@isnumeric, cols(isalgo)))
    error('threshold_sweep:nonNumeric', ['At least one algo_ column is ' ...
        'non-numeric. Check the CSV for "NA", blank rows, or Excel error ' ...
        'codes such as #REF! or #DIV/0!.']);
end
Y = [cols{isalgo}];
labels = regexprep(algoNames, '^algo_', '', 'ignorecase');
% The Source column, when present, names each instance's benchmark class or
% generator. It is the population the degeneracy guard checks against -- a
% threshold that merely re-derives it has sorted instances, not algorithms.
issource = strcmpi(names, 'source');
groups = {};
if any(issource)
    [~, ~, cells] = readCsv(csvPath);
    groups = strtrim(cells(:, find(issource, 1)));
end
end

function v = extremum(x, fn)
% extremum  min or max of x, or NaN when x is empty.
if isempty(x)
    v = NaN;
else
    v = fn(x);
end
end

function [gidx, ngroups, ordered] = groupIndex(groups)
% groupIndex  Map group labels to an ordered index.
%   Benchmark classes usually carry a trailing number whose numeric order is
%   the meaningful one -- alphabetical would put set10 before set2.
%   ordered is false when the labels are not numbered, or when two families
%   share a numbering (setA02 and setB02 both parse to 2); groupRho is then
%   meaningless and returned as NaN, leaving groupSpread as the guard.
if isempty(groups)
    gidx = [];
    ngroups = 0;
    ordered = false;
    return;
end
[uniq, ~, idx] = unique(groups(:));
ngroups = numel(uniq);
nums = nan(ngroups, 1);
for k = 1:ngroups
    tok = regexp(uniq{k}, '(\d+)\s*$', 'tokens', 'once');
    if ~isempty(tok)
        nums(k) = str2double(tok{1});
    end
end
if all(~isnan(nums)) && numel(unique(nums)) == ngroups
    [~, order] = sort(nums);
    remap = zeros(ngroups, 1);
    remap(order) = 1:ngroups;
    gidx = remap(idx);
    ordered = true;
else
    gidx = idx;
    ordered = false;
end
end

function rho = spearmanRho(x, y)
% spearmanRho  Spearman rank correlation, average ranks for ties.
if numel(x) < 3
    rho = NaN;
    return;
end
c = corrcoef(averageRanks(x(:)), averageRanks(y(:)));
rho = c(1, 2);
if ~isfinite(rho)
    rho = NaN;
end
end

function rk = averageRanks(x)
% averageRanks  1-based ranks, tied values sharing their average rank.
[s, order] = sort(x);
rk = zeros(size(x));
n = numel(x);
i = 1;
while i <= n
    j = i;
    while j < n && s(j + 1) == s(i)
        j = j + 1;
    end
    rk(order(i:j)) = (i + j) / 2;
    i = j + 1;
end
end

function v = meanPairwiseJaccard(lab)
% meanPairwiseJaccard  Mean overlap between algorithms' good-sets.
%   Lower means the portfolio is more complementary at this threshold. Pairs
%   whose union is empty carry no information and are skipped.
n = size(lab, 2);
if n < 2
    v = NaN;
    return;
end
vals = [];
for a = 1:n - 1
    for b = a + 1:n
        u = nnz(lab(:, a) | lab(:, b));
        if u
            vals(end + 1) = nnz(lab(:, a) & lab(:, b)) / u; %#ok<AGROW>
        end
    end
end
if isempty(vals)
    v = NaN;
else
    v = mean(vals);
end
end

function [slope, width] = plateauStats(grid, gRate, e, gAtE, tol, rateLow, rateHigh)
% plateauStats  Local slope and plateau width around a candidate rung.
%   Width is measured on the dense curve: the widest contiguous interval
%   containing e over which the good-rate stays within tol of its value at
%   e. A wide plateau is the evidence that the exact value is not doing the
%   work.
slope = NaN;
width = NaN;
if e < grid(1) || e > grid(end)
    return;
end
[~, i] = min(abs(grid - e));
if i > 1 && i < numel(grid)
    slope = (gRate(i + 1) - gRate(i - 1)) / (grid(i + 1) - grid(i - 1));
end
% A flat curve at a saturated good-rate is degeneracy, not robustness: the
% labels have stopped changing, so of course epsilon can move freely. Report
% no width rather than one that would be misread as evidence of stability.
if gAtE < rateLow || gAtE > rateHigh
    return;
end
ok = abs(gRate - gAtE) <= tol;
if ~ok(i)
    width = 0;
    return;
end
lo = i;
hi = i;
while lo > 1 && ok(lo - 1)
    lo = lo - 1;
end
while hi < numel(grid) && ok(hi + 1)
    hi = hi + 1;
end
width = grid(hi) - grid(lo);
end

function e = ladder125(lo, hi)
% ladder125  Fallback 1-2-5 sweep grid spanning the data range.
%   Exists so the function runs without a caller-supplied ladder. It is a
%   sweep convenience, not a proposal -- a real ladder comes from the
%   metric's units and the domain's conventions (see the skill's "Generate
%   the candidate ladder", section 3).
if ~isfinite(lo) || ~isfinite(hi) || hi <= 0
    e = [0 0.01 0.02 0.05 0.1 0.2 0.5 1];
    return;
end
lo = max(lo, hi / 1e4);
decades = floor(log10(lo)):ceil(log10(hi));
e = unique(reshape([1; 2; 5] * 10 .^ decades, 1, []));
e = e(e >= lo / 2 & e <= hi * 1.5);
if isempty(e)
    e = linspace(lo, hi, 7);
end
end

% =========================================================================
% Report
% =========================================================================
function printReport(r)
% printReport  Print the sweep in the order the skill reads it.
fprintf('== data (threshold-independent; read this first) ==\n');
keys = fieldnames(r.data);
for k = 1:numel(keys)
    key = keys{k};
    v = r.data.(key);
    if strcmp(key, 'algoLabels')
        shown = strjoin(v, ', ');
    elseif strcmp(key, 'safeRungInterval')
        if isempty(v)
            continue;
        elseif ischar(v)
            fprintf('  %-24spoint mass sits on the GOOD side of any value: if it is a sentinel, no epsilon labels it bad\n', key);
            continue;
        elseif ~isfinite(v(1)) || ~isfinite(v(2))
            shown = 'no genuine values on the safe side of the sentinel';
        else
            shown = sprintf(['a threshold strictly between %g and %g labels the ' ...
                'sentinel bad and every genuine run good'], v(1), v(2));
        end
    elseif islogical(v)
        shown = mat2str(v);
    else
        shown = sprintf('%.10g', v);
    end
    fprintf('  %-24s%s\n', key, shown);
end
fprintf('\n== per algorithm ==\n');
printTable(r.perAlgo, @(x) sprintf('%.4g', x));
% Two blocks rather than one, both keyed by epsilon so rows still line up
% across them. Every r.sweep column the skill's references name as a
% diagnostic has to be visible here, and all of them on one line is 200
% characters. epsilon gets %g, not the fixed 3 decimals the diagnostics get:
% candidate values span 0.0001 to 3600 depending on the metric and the mode,
% and rounding would print 0 and 0.0001 identically.
fmt = @(x) sprintf('%.3f', x);
fprintf('\n== sweep: labels ==\n');
printTable(r.sweep(:, {'epsilon', 'goodRate', 'algoRateMin', 'algoRateMax', ...
    'fracAllGood', 'fracAllBad', 'nDegenerateAlgos', 'nNearDegenerate'}), fmt);
fprintf('\n== sweep: selection and shape ==\n');
printTable(r.sweep(:, {'epsilon', 'meanJaccard', 'uniqueBestFrac', 'nPassengers', ...
    'dGoodRate', 'plateauWidth', 'groupSpread', 'groupRho', 'screenScore'}), fmt);
fprintf('\n== unique contribution per algorithm ==\n');
for ri = 1:height(r.sweep)
    parts = cell(1, numel(r.data.algoLabels));
    for a = 1:numel(parts)
        parts{a} = sprintf('%s=%.3f', r.data.algoLabels{a}, r.uniqueByAlgo{ri}(a));
    end
    fprintf('  eps=%-12g%s\n', r.sweep.epsilon(ri), strjoin(parts, '  '));
end
end

function printTable(T, fmt)
% printTable  Right-aligned columns; epsilon as %g, counts as integers.
names = T.Properties.VariableNames;
nrow = height(T);
cells = cell(nrow + 1, numel(names));
for c = 1:numel(names)
    cells{1, c} = names{c};
    col = T.(names{c});
    for i = 1:nrow
        if iscell(col)
            s = col{i};
        elseif strcmp(names{c}, 'epsilon')
            s = sprintf('%g', col(i));
        elseif isnan(col(i))
            s = 'NaN';
        elseif any(strcmp(names{c}, {'nDegenerateAlgos', 'nNearDegenerate', 'nPassengers'}))
            s = sprintf('%d', col(i));
        else
            s = fmt(col(i));
        end
        cells{i + 1, c} = s;
    end
end
widths = max(cellfun(@numel, cells), [], 1);
for i = 1:nrow + 1
    line = '';
    for c = 1:numel(names)
        line = [line, sprintf('%*s', widths(c) + 2 * (c > 1), cells{i, c})]; %#ok<AGROW>
    end
    fprintf('%s\n', line);
end
end

function [names, cols, cells] = readCsv(path)
% readCsv  Every column of a CSV, numeric where each cell is a number or blank.
%   Read by hand rather than with readtable, which guesses the delimiter
%   (MATLAB picks '_' on some narrow files), renames headers that are not
%   valid identifiers, and in Octave fails on any empty cell. The delimiter
%   is whichever of comma, semicolon or tab the header holds most of. A
%   blank cell, or NaN spelled out, is read as NaN; any other cell that is
%   not a number makes its column text. cells keeps every cell as written.
if exist(path, 'file') ~= 2
    error('readCsv:cannotOpen', 'cannot open %s', path);
end
txt = fileread(path);
if ~isempty(txt) && double(txt(1)) == 65279, txt = txt(2:end); end  % BOM
lines = regexp(txt, '\r?\n', 'split');
lines = lines(~cellfun(@(l) isempty(strtrim(l)), lines));
if isempty(lines)
    error('readCsv:empty', '%s is empty', path);
end
delims = {',', ';', sprintf('\t')};
counts = cellfun(@(d) sum(lines{1} == d), delims);
[top, best] = max(counts);
delim = ',';
if top > 0, delim = delims{best}; end
names = splitLine(lines{1}, delim);
cells = cell(numel(lines) - 1, numel(names));
for i = 2:numel(lines)
    f = splitLine(lines{i}, delim);
    f(end + 1:numel(names)) = {''};
    cells(i - 1, :) = f(1:numel(names));
end
cols = cell(1, numel(names));
for k = 1:numel(names)
    s = strtrim(cells(:, k));
    x = str2double(s);
    blank = cellfun(@isempty, s);
    lit = ~cellfun(@isempty, regexpi(s, '^[+-]?(nan|inf|infinity)$', 'once'));
    if all((~isnan(x) & imag(x) == 0) | blank | lit)
        x = real(x);
        x(blank) = NaN;
        cols{k} = x;
    else
        cols{k} = cells(:, k);
    end
end
end

function f = splitLine(line, delim)
% splitLine  Split one CSV line, honouring double-quoted fields.
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
