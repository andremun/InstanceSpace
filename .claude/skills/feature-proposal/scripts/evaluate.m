function varargout = evaluate(featuresCsv, performanceCsv, varargin)
% evaluate  Run the feature-proposal checks over a set of candidate features.
%
%   evaluate('proposed.csv', 'perf.csv')
%   evaluate('proposed.csv', 'perf.csv', 'Existing', 'current.csv', ...
%            'MaxPerf', false, 'AbsPerf', false, 'Epsilon', 0.05)
%   r = evaluate(...)          % returns the results instead of printing them
%
%   CSV shape: the first column is the instance id; any further non-numeric
%   column is ignored; every numeric column is a feature (or, in the
%   performance file, an algorithm). Ids must match across files.
%
%   The correlation check reproduces SIFTED's own filter, which is what
%   decides whether a candidate survives selection. So the features and
%   performance go through the toolkit's own PRELIM, the correlation is the
%   same corr(..., 'rows', 'pairwise') call SIFTED makes, and the clusters
%   come from kmeans with SIFTED's settings and seed. Every default is read
%   from ISAdefaults. Against a built space, pass that space's
%   model.opts.perf and model.opts.sifted values instead of trusting the
%   defaults -- the run may have set them differently.
%
%   To audit a feature set rather than screen candidates against one, pass
%   the set itself as the features file and leave Existing off. The
%   degeneracy, correlation and clustering reports then describe that set.
%
%   Options (name-value pairs; defaults from ISAdefaults):
%     Existing    CSV of the feature set already in use (enables the overlap,
%                 incremental-value and joint clustering reports)
%     MaxPerf     opts.perf.MaxPerf: true = higher is better
%     AbsPerf     opts.perf.AbsPerf: true = fixed threshold, false = relative
%     Epsilon     opts.perf.epsilon: the good-performance threshold
%     Rho         opts.sifted.rho: minimum |correlation| to keep
%     Pval        opts.sifted.pval: significance level for that correlation
%     K           opts.sifted.K: how many clusters SIFTED forms
%     Print       print the report (default: true when no output is requested)
%
%   Runs in MATLAB and in GNU Octave with the statistics and datatypes
%   packages loaded, as the toolkit itself needs. Run it with the toolkit on
%   the path: it calls PRELIM, ISAdefaults and isacompat.
%
%   Companion function for the feature-proposal skill.
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

MIN_ROWS = 3;          % fewest observed rows a column needs to be correlated
MIN_DISTINCT = 2;      % fewest distinct values a column or label needs to correlate
NEAR_CONSTANT = 0.90;  % share of rows at the modal value that earns a WARN
CV_FOLDS = 5;          % folds for the incremental-value check

d = ISAdefaults(struct());
p = inputParser;
p.FunctionName = 'evaluate';
addRequired(p, 'featuresCsv', @(x) ischar(x) || isstring(x));
addRequired(p, 'performanceCsv', @(x) ischar(x) || isstring(x));
addParameter(p, 'Existing', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'MaxPerf', d.perf.MaxPerf, @isFlag);
addParameter(p, 'AbsPerf', d.perf.AbsPerf, @isFlag);
addParameter(p, 'Epsilon', d.perf.epsilon, @(x) isscalar(x) && isnumeric(x));
addParameter(p, 'Rho', d.sifted.rho, @(x) isscalar(x) && isnumeric(x));
addParameter(p, 'Pval', d.sifted.pval, @(x) isscalar(x) && isnumeric(x));
addParameter(p, 'K', d.sifted.K, @(x) isscalar(x) && isnumeric(x) && x >= 1);
addParameter(p, 'Print', nargout == 0, @isFlag);
parse(p, featuresCsv, performanceCsv, varargin{:});
a = p.Results;
existing = char(a.Existing);
say = makePrinter(logical(a.Print));

[ids, names, featRaw] = readNumericCsv(char(a.featuresCsv));
[pids, algos, perfRaw] = readNumericCsv(char(a.performanceCsv));
perfRaw = alignRows(ids, pids, perfRaw, char(a.performanceCsv));
[n, m] = size(featRaw);
nAlgos = numel(algos);

% ---- the good/bad rule, and the target SIFTED would correlate against ----
popts = struct('MaxPerf', logical(a.MaxPerf), 'AbsPerf', logical(a.AbsPerf), ...
    'epsilon', a.Epsilon, 'betaThreshold', d.perf.betaThreshold, 'auto', false, ...
    'bound', false, 'norm', false, 'iqrMultiplier', d.prelim.iqrMultiplier);
yBin = runPrelim(ones(n, 1), perfRaw, popts);
yBin = double(yBin);
if popts.AbsPerf
    kind = 'absolute';
    target = 'raw performance';
else
    kind = 'relative';
    target = 'deficit from the best';
end
say('n = %d instances, %d candidate features, %d algorithms\n', n, m, nAlgos);
say('good rule: %s, epsilon %g, MaxPerf %s; the correlation target is the %s\n', ...
    kind, a.Epsilon, mat2str(popts.MaxPerf), target);
rates = cell(1, nAlgos);
for k = 1:nAlgos
    rates{k} = sprintf('%s %.0f%%', algos{k}, 100 * mean(yBin(:, k)));
end
say('  good rate: %s\n', strjoin(rates, ', '));
dead = mean(sum(yBin, 2) == 0);
note = '';
if ~popts.AbsPerf
    note = '  <- these are ties in the deficit target';
end
say('  instances no algorithm is good on: %.0f%%%s\n\n', 100 * dead, note);

% The build removes an algorithm that is good on no instance before PRELIM
% normalises anything (InstanceSpace.m:636-645), so it is not in the matrix
% SIFTED correlates against. Leave it out here too.
never = sum(yBin, 1) == 0;
if all(never)
    error('evaluate:noPortfolio', ['no algorithm is good on any instance under ' ...
        'this rule, so the build would have no portfolio left']);
end
if any(never)
    say('  never good, so the build removes them before anything else runs: %s. They are left out below.\n\n', ...
        strjoin(algos(never), ', '));
    perfRaw = perfRaw(:, ~never);
    yBin = yBin(:, ~never);
    algos = algos(~never);
    nAlgos = numel(algos);
end

% ---- reject 1: degeneracy -----------------------------------------------
say('== reject 1: degeneracy ==\n');
say('%-28s %8s %6s %8s %9s\n', 'feature', 'distinct', 'mode', 'missing', 'minority');
alive = true(1, m);
[~, spread] = isacompat.columnQuartiles(featRaw);
degen = cell(m, 1);
for j = 1:m
    v = featRaw(:, j);
    obs = v(~isnan(v));
    [vals, ~, ic] = unique(obs);
    counts = accumarray(ic(:), 1);
    k = numel(vals);
    if isempty(obs)
        mode = 1;
        top = 0;
    else
        top = max(counts);
        mode = top / numel(obs);
    end
    missing = mean(isnan(v));
    note = '';
    minority = '        -';
    if k <= 1
        alive(j) = false;
        note = '  FAIL constant';
    elseif missing >= d.prelim.nanThreshold
        alive(j) = false;
        note = sprintf('  FAIL missing >= %.0f%%, INIT would drop it', 100 * d.prelim.nanThreshold);
    elseif spread(j) == 0
        alive(j) = false;
        minority = sprintf('%9d', numel(obs) - top);
        note = '  FAIL zero IQR, PRELIM''s bounding collapses it to a constant';
    elseif mode >= NEAR_CONSTANT
        minority = sprintf('%9d', numel(obs) - top);
        note = sprintf('  WARN near-constant, %d rows carry the signal', numel(obs) - top);
    end
    degen{j} = strtrim(note);
    say('%-28s %8d %6.2f %8.2f %s%s\n', lab(names{j}, 28), k, mode, missing, minority, note);
end
say(['  %d/%d pass. A WARN passes through to the next check; all it costs is\n' ...
     '  sample size, since the correlation rests on the minority rows alone.\n'], sum(alive), m);
say(['  A zero IQR is different and is fatal as the pipeline is configured: the middle 50%%\n' ...
     '  of the column is a single value, so median +/- %g*IQR has zero width and PRELIM''s\n' ...
     '  bounding step (core/PRELIM.m:225) clips every row to that value. SIFTED then sees a\n' ...
     '  constant column, scores it 0 (core/SIFTED.m:117) and drops it. Any binary column\n' ...
     '  whose minority is under a quarter of the rows lands here. To keep such a feature,\n' ...
     '  set opts.bound.flag false, or redefine it so the middle 50%% is not one value.\n\n'], ...
     d.prelim.iqrMultiplier);

% ---- reject 2: SIFTED's correlation filter ------------------------------
ppopts = popts;
ppopts.auto = d.auto.preproc;
ppopts.bound = d.bound.flag;
ppopts.norm = d.norm.flag;
[feat, yProc] = runPrelim(featRaw, perfRaw, ppopts, 'processed');
[rMat, pMat] = corr(feat, yProc, 'rows', 'pairwise');
rMat(isnan(rMat)) = 0;
pMat(isnan(pMat)) = 1;

say('== reject 2: Pearson on processed values, |rho| >= %g at p <= %g ==\n', a.Rho, a.Pval);
say('this is core/SIFTED.m:115 -- same target, same transform, same thresholds, no correction\n');
say('%-28s %s %7s\n', 'feature', strjoin(cellfun(@(c) sprintf('%10s', clip(c, 9)), algos, 'UniformOutput', false), ' '), 'spread');
hitMat = abs(rMat) >= a.Rho & pMat <= a.Pval;
passes = alive & any(hitMat, 2)';
for j = 1:m
    cells = cell(1, nAlgos);
    for k = 1:nAlgos
        star = ' ';
        if hitMat(j, k), star = '*'; end
        cells{k} = sprintf('%+9.3f%s', rMat(j, k), star);
    end
    drop = '';
    if ~passes(j), drop = '  drop'; end
    say('%-28s %s %7.3f%s\n', lab(names{j}, 28), strjoin(cells, ' '), max(rMat(j, :)) - min(rMat(j, :)), drop);
end
say(['  %d/%d pass. ''*'' marks an algorithm this feature clears the floor on;\n' ...
     '  spread is the range across the portfolio, high meaning the feature separates\n' ...
     '  algorithms rather than only tracking difficulty.\n'], sum(passes), m);
say(['  The floor is sufficient, not necessary: SIFTED also keeps whichever feature is\n' ...
     '  most correlated with each algorithm whatever its value (core/SIFTED.m:121), so in a\n' ...
     '  real run a candidate below the floor can still survive by being the best there is.\n']);

curved = zeros(0, 3);
for j = 1:m
    if passes(j) || ~alive(j), continue; end
    for k = 1:nAlgos
        if hitMat(j, k), continue; end
        [rq, pq] = quadFit(feat(:, j), yProc(:, k));
        if rq >= a.Rho && pq <= a.Pval && rq > abs(rMat(j, k))
            curved(end + 1, :) = [j, k, rq]; %#ok<AGROW>
        end
    end
end
if ~isempty(curved)
    say(['\n  NON-MONOTONE SIGNAL in candidates the check above drops. A quadratic reaches\n' ...
         '  the same floor where Pearson does not, so the relationship is real and runs\n' ...
         '  both ways around a turning point:\n']);
    [~, order] = sort(-curved(:, 3));
    for t = order'
        j = curved(t, 1); k = curved(t, 2);
        say('    %-28s %-20s Pearson %+.3f, quadratic R %.3f\n', lab(names{j}, 28), ...
            lab(algos{k}, 20), rMat(j, k), curved(t, 3));
    end
    say(['  The pipeline is Pearson throughout, so it will not find these as they\n' ...
         '  stand. Redefine the feature so the relationship runs one way -- the distance\n' ...
         '  from the turning point rather than the raw value -- and run it again.\n']);
end
say('\n');
if ~any(passes)
    error('evaluate:nothingSurvives', 'nothing survives the correlation check');
end
survivors = find(passes);

% ---- rank ---------------------------------------------------------------
rbMat = zeros(m, nAlgos);
for j = 1:m
    ok = ~isnan(featRaw(:, j));
    for k = 1:nAlgos
        if sum(ok) < MIN_ROWS || numel(unique(yBin(ok, k))) < MIN_DISTINCT || ...
                numel(unique(featRaw(ok, j))) < MIN_DISTINCT
            continue;
        end
        c = corrcoef(yBin(ok, k), featRaw(ok, j));
        if isfinite(c(1, 2)), rbMat(j, k) = c(1, 2); end
    end
end

existProc = [];
enames = {};
overlap = table();
if ~isempty(existing)
    [eids, enames, existRaw] = readNumericCsv(existing);
    existRaw = alignRows(ids, eids, existRaw, existing);
    existProc = runPrelim(existRaw, perfRaw, ppopts, 'processed');
    near = cell(numel(survivors), 1);
    rNear = zeros(numel(survivors), 1);
    for s = 1:numel(survivors)
        rr = abs(corr(feat(:, survivors(s)), existProc, 'rows', 'pairwise'));
        rr(isnan(rr)) = 0;
        [rNear(s), kNear] = max(rr);
        near{s} = enames{kNear};
    end
    overlap = table(names(survivors)', rNear, near, 'VariableNames', {'feature', 'maxAbsR', 'nearestExisting'});
end

say('== rank: point-biserial against the binary label, per algorithm ==\n');
say('what PYTHIA predicts and TRACE draws footprints from, so it can order the\n');
say('survivors differently from the check above\n');
say('%-28s %s %7s\n', 'feature', strjoin(cellfun(@(c) sprintf('%10s', clip(c, 9)), algos, 'UniformOutput', false), ' '), 'spread');
[~, order] = sort(-max(abs(rbMat(survivors, :)), [], 2));
for j = survivors(order)
    say('%-28s %s %7.3f\n', lab(names{j}, 28), strjoin(arrayfun(@(x) sprintf('%+10.3f', x), rbMat(j, :), 'UniformOutput', false), ' '), ...
        max(rbMat(j, :)) - min(rbMat(j, :)));
end
say('\n');

if ~isempty(existing)
    say('== rank: overlap with %d existing features ==\n', numel(enames));
    say(['whether the set already measures this. A duplicate does no harm, since SIFTED\n' ...
         'clusters the two and keeps one, but it buys nothing either\n']);
    say('%-28s %7s  %-28s\n', 'feature', 'max|r|', 'nearest existing');
    [~, order] = sort(-overlap.maxAbsR);
    for s = order'
        say('%-28s %7.3f  %-28s\n', lab(overlap.feature{s}, 28), overlap.maxAbsR(s), lab(overlap.nearestExisting{s}, 28));
    end
    say('\n');
end

poolNames = names(survivors);
poolKind = repmat({'candidate'}, 1, numel(survivors));
poolMat = feat(:, survivors);
nLeftOut = 0;
if ~isempty(existing)
    % SIFTED clusters only the features that pass its correlation filter, so an
    % existing feature that would fail it -- a constant after PRELIM's
    % bounding, say -- never reaches the clustering. Leave it out here too.
    [rE, pE] = corr(existProc, yProc, 'rows', 'pairwise');
    rE(isnan(rE)) = 0;
    pE(isnan(pE)) = 1;
    eKeep = any(abs(rE) >= a.Rho & pE <= a.Pval, 2)';
    nLeftOut = sum(~eKeep);
    poolNames = [poolNames, enames(eKeep)];
    poolKind = [poolKind, repmat({'existing'}, 1, sum(eKeep))];
    poolMat = [poolMat, existProc(:, eKeep)];
end
clusters = {};
say('== rank: the clusters SIFTED would form ==\n');
if nLeftOut > 0
    say('  %d existing feature(s) left out: they would not pass SIFTED''s correlation filter.\n', nLeftOut);
end
npool = numel(poolNames);
if npool <= max(a.K, 3)
    plural = 's';
    if npool == 1, plural = ''; end
    say(['  %d column%s is not more than K = %d, so SIFTED skips the clustering\n' ...
         '  and keeps all of them (core/SIFTED.m:140).\n\n'], npool, plural, a.K);
else
    grp = siftedClusters(poolMat, a.K, d.sifted);
    say(['  k-means on correlation distance into opts.sifted.K clusters, with SIFTED''s\n' ...
         '  own settings and seed (core/SIFTED.m:176). It then keeps exactly one column from\n' ...
         '  each (core/SIFTED.m:200-204), so columns sharing a cluster are competing and\n' ...
         '  all but one of them will be dropped whatever you add to the set.\n']);
    for g = 1:a.K
        members = find(grp == g);
        if isempty(members), continue; end
        tags = cell(1, numel(members));
        for t = 1:numel(members)
            tags{t} = lab(poolNames{members(t)}, 24);
            if ~isempty(existing) && strcmp(poolKind{members(t)}, 'existing')
                tags{t} = [tags{t} ' (existing)'];
            end
        end
        clusters{end + 1} = poolNames(members); %#ok<AGROW>
        say('    cluster %d: %s\n', g, strjoin(tags, ', '));
    end
    nused = numel(unique(grp(~isnan(grp))));
    say('  %d columns in %d clusters, so the set comes out of selection\n  at %d columns however many go in.\n\n', ...
        npool, nused, nused);
end

incremental = table();
if ~isempty(existing)
    xExisting = existProc;
    xCombined = [existProc, feat(:, survivors)];
    say('== rank: incremental value, %d-fold CV AUC, existing vs existing + survivors ==\n', CV_FOLDS);
    say('%-16s %9s %11s %8s %9s\n', 'algorithm', 'existing', '+survivors', 'gain', 'headroom');
    rows = {};
    for k = 1:nAlgos
        y = yBin(:, k);
        nMinority = min(sum(y), numel(y) - sum(y));
        if numel(unique(y)) < MIN_DISTINCT
            say('%-16s  single class, skipped\n', lab(algos{k}, 16));
            continue;
        end
        if nMinority < CV_FOLDS
            say('%-16s  %d in the smaller class, too few for %d-fold, skipped\n', lab(algos{k}, 16), nMinority, CV_FOLDS);
            continue;
        end
        s1 = cvAuc(xExisting, y, CV_FOLDS);
        s2 = cvAuc(xCombined, y, CV_FOLDS);
        rows(end + 1, :) = {algos{k}, s1, s2}; %#ok<AGROW>
        say('%-16s %9.3f %11.3f %+8.3f %9.3f\n', lab(algos{k}, 16), s1, s2, s2 - s1, 1 - s1);
    end
    say(['  read the gain against the headroom beside it: a small gain on an already-high\n' ...
         '  baseline is a ceiling effect, not evidence the candidates add nothing\n\n']);
    if ~isempty(rows)
        incremental = cell2table(rows, 'VariableNames', {'algorithm', 'aucExisting', 'aucWithSurvivors'});
    end
end

say('== rank: ordered ==\n');
say('on |rho| against the performance values, because that is what decides whether\n');
say('SIFTED keeps the feature. ''on'' names the algorithm each figure belongs to; the\n');
say('two can be different algorithms, and either column can order the set its own\n');
say('way. Discount a candidate whose overlap above is near 1 however high its rho --\n');
say('the order below does not apply that discount for you.\n');
say('%-28s %8s %-11s %9s %-11s %7s\n', 'feature', '|rho| Y', ' on', '|r| Ybin', ' on', 'spread');
[~, order] = sort(-max(abs(rMat(survivors, :)), [], 2));
for j = survivors(order)
    [ry, ky] = max(abs(rMat(j, :)));
    [rb, kb] = max(abs(rbMat(j, :)));
    say('%-28s %8.3f %-11s %9.3f %-11s %7.3f\n', lab(names{j}, 28), ry, [' ' lab(algos{ky}, 10)], ...
        rb, [' ' lab(algos{kb}, 10)], max(rMat(j, :)) - min(rMat(j, :)));
end

% ---- yield --------------------------------------------------------------
say('\n== yield ==\n');
say('  %d of %d candidates evaluated in this run survive both checks.\n', numel(survivors), m);
say('  Report that against every candidate ever evaluated, not just this run''s.\n');

if nargout > 0
    r = struct();
    r.features = names;
    r.algorithms = algos;
    r.neverGood = never;
    r.degeneracy = table(names(:), alive(:), degen, 'VariableNames', {'feature', 'alive', 'note'});
    r.rho = rMat;
    r.pval = pMat;
    r.passes = passes;
    r.survivors = names(survivors);
    r.pointBiserial = rbMat;
    r.overlap = overlap;
    r.clusters = clusters;
    r.incremental = incremental;
    varargout{1} = r;
end
end

% =========================================================================
% Helpers
% =========================================================================
function tf = isFlag(x)
tf = isscalar(x) && (islogical(x) || isnumeric(x));
end

function say = makePrinter(enabled)
% makePrinter  fprintf when printing is on, a no-op otherwise.
if enabled
    say = @(varargin) fprintf(varargin{:});
else
    say = @(varargin) [];
end
end

function s = lab(name, w)
% lab  Feature names run long; keep the tail, which is what distinguishes them.
if numel(name) <= w
    s = name;
else
    s = ['..' name(end - w + 3:end)];
end
end

function s = clip(name, w)
% clip  The first w characters of a name.
s = name(1:min(end, w));
end

function [ids, names, X] = readNumericCsv(path)
% readNumericCsv  Id column, then every numeric column as a matrix.
[header, cols, cells] = readCsv(path);
ids = cells(:, 1);
keep = false(1, numel(header));
for c = 2:numel(header)
    keep(c) = isnumeric(cols{c});
end
if ~any(keep)
    error('evaluate:noColumns', '%s: no numeric columns after the id column', path);
end
names = header(keep);
X = [cols{keep}];
end

function Xb = alignRows(idsA, idsB, Xb, what)
% alignRows  Reorder the rows of Xb, keyed by idsB, into the order of idsA.
[found, pos] = ismember(idsA, idsB);
if ~all(found)
    missing = idsA(~found);
    error('evaluate:idsAbsent', '%s: %d ids absent, e.g. %s', what, numel(missing), ...
        strjoin(missing(1:min(3, end))', ', '));
end
Xb = Xb(pos, :);
end

function varargout = runPrelim(X, Y, opts, mode)
% runPrelim  Call the toolkit's PRELIM without its progress messages.
%   Returns Ybin, or with mode 'processed' the processed X and Y.
if nargin < 4, mode = 'labels'; end
[~, Xp, Yp, out] = evalcPrelim(X, Y, opts);
if strcmp(mode, 'processed')
    varargout = {Xp, Yp};
else
    varargout = {out.Ybin};
end
end

function [txt, Xp, Yp, out] = evalcPrelim(X, Y, opts) %#ok<INUSD,STOUT> used inside evalc
txt = evalc('[Xp, Yp, out] = PRELIM(X, Y, opts);');
end

function [r, p] = quadFit(x, y)
% quadFit  Multiple R for y ~ 1 + x + x^2 and its F-test p, pairwise-complete.
%   Pearson scores a U-shape near zero, so a feature that predicts
%   performance on both sides of a turning point fails the correlation check
%   at nothing. This is what finds those before they are discarded in silence.
r = 0;
p = 1;
ok = ~isnan(x) & ~isnan(y);
if sum(ok) < 10 || numel(unique(x(ok))) < 3, return; end
xs = x(ok);
ys = y(ok);
if all(abs(ys - ys(1)) <= 1e-8 + 1e-5 * abs(ys(1))), return; end
design = [ones(numel(xs), 1), xs, xs .^ 2];
beta = design \ ys;
sstot = sum((ys - mean(ys)) .^ 2);
if sstot <= 0, return; end
r2 = max(0, 1 - sum((ys - design * beta) .^ 2) / sstot);
dfres = numel(xs) - 3;
if dfres <= 0, return; end
if r2 >= 1
    r = 1;
    p = 0;
    return;
end
f = (r2 / 2) / ((1 - r2) / dfres);
r = sqrt(r2);
p = betainc(dfres / (dfres + 2 * f), dfres / 2, 1);  % upper tail of F(2, dfres)
end

function grp = siftedClusters(M, K, s)
% siftedClusters  Group columns of M the way SIFTED does (core/SIFTED.m:176).
%   kmeans on correlation distance, over features as points, with SIFTED's
%   MaxIter, Replicates and seed. A missing value is replaced by its feature's
%   mean first, since kmeans would otherwise drop the whole feature.
M = M';
mu = mean(M, 2, 'omitnan');
[ri, ci] = find(isnan(M));
M(sub2ind(size(M), ri, ci)) = mu(ri);
prev = rng;
restore = onCleanup(@() rng(prev));
rng(s.seed, 'twister');
args = {'Distance', 'correlation', 'MaxIter', s.MaxIter, 'Replicates', s.Replicates};
if ~isacompat.isOctave()
    args = [args, {'OnlinePhase', 'on'}];
end
grp = kmeans(M, K, args{:});
end

function auc = cvAuc(X, y, folds)
% cvAuc  Mean K-fold ROC AUC of an L2 logistic regression on X.
%   Folds are stratified and deterministic: each class's instances are dealt
%   round-robin in file order, so MATLAB and Octave give the same result.
%   Missing values take the training fold's mean.
fold = zeros(numel(y), 1);
for cls = [0 1]
    idx = find(y == cls);
    fold(idx) = mod(0:numel(idx) - 1, folds) + 1;
end
aucs = zeros(1, folds);
for f = 1:folds
    tr = fold ~= f;
    te = ~tr;
    mu = mean(X(tr, :), 1, 'omitnan');
    mu(isnan(mu)) = 0;
    Xtr = fillMissing(X(tr, :), mu);
    Xte = fillMissing(X(te, :), mu);
    w = fitLogistic(Xtr, y(tr));
    aucs(f) = rocAuc([ones(sum(te), 1), Xte] * w, y(te));
end
auc = mean(aucs);
end

function X = fillMissing(X, mu)
[ri, ci] = find(isnan(X));
X(sub2ind(size(X), ri, ci)) = mu(ci);
end

function w = fitLogistic(X, y)
% fitLogistic  Newton's method for logistic loss + 0.5*||w||^2, intercept unpenalised.
A = [ones(size(X, 1), 1), X];
w = zeros(size(A, 2), 1);
pen = eye(size(A, 2));
pen(1, 1) = 0;
for it = 1:100
    pr = 1 ./ (1 + exp(-A * w));
    g = A' * (pr - y) + pen * w;
    H = A' * bsxfun(@times, A, pr .* (1 - pr)) + pen + 1e-10 * eye(size(A, 2));
    step = H \ g;
    w = w - step;
    if max(abs(step)) < 1e-10, break; end
end
end

function auc = rocAuc(score, y)
% rocAuc  Area under the ROC curve via the rank-sum statistic, ties averaged.
pos = y == 1;
np = sum(pos);
nn = sum(~pos);
if np == 0 || nn == 0
    auc = NaN;
    return;
end
[s, order] = sort(score);
rk = zeros(size(score));
i = 1;
n = numel(score);
while i <= n
    j = i;
    while j < n && s(j + 1) == s(i)
        j = j + 1;
    end
    rk(order(i:j)) = (i + j) / 2;
    i = j + 1;
end
auc = (sum(rk(pos)) - np * (np + 1) / 2) / (np * nn);
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
