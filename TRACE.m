function out = TRACE(Z, Ybin, Yhat, P, beta, algolabels, opts, trainedTrace)
% TRACE  Compute or evaluate algorithm footprints in the instance space.
%
%   Training mode (7 args):
%     out = TRACE(Z, Ybin, Yhat, P, beta, algolabels, opts)
%     Builds footprints using TRACE3 (default) or the legacy DBSCAN algorithm.
%
%   Evaluation mode (8 args):
%     out = TRACE(Z, Ybin, Yhat, P, beta, algolabels, opts, trainedTrace)
%     Re-evaluates trained footprints on new instances; polygons are not rebuilt.
%
%   Inputs
%     Z          - (ninst x ndim) projected instance coordinates (2D or 3D)
%     Ybin       - (ninst x nalgos) true binary performance labels
%     Yhat       - (ninst x nalgos) PYTHIA predicted labels; [] if unavailable
%     P          - (ninst x 1) best-algorithm index per instance
%     beta       - (ninst x 1) logical: instance is easy (has a good algorithm)
%     algolabels - (1 x nalgos) cell array of algorithm name strings
%     opts       - struct with fields (see ISAdefaults for defaults):
%                    method        'trace3' (default) or 'legacy'
%                    PI            minimum purity threshold (0.6)
%                    nn            (unused; reserved for future KNN fallback)
%                    prior         (unused; reserved for future KNN fallback)
%                  Note: in TRACE3, if Yhat is empty, Zu = {yi=1} directly;
%                  there is no KNN fallback in this version.
%                    minInstances  minimum instances for a valid footprint (4)
%                    minAreaFrac   footprint must exceed this fraction of space (0.01)
%                    contra        contradiction removal; legacy only; default true
%                                  when method='legacy', false otherwise
%     trainedTrace - (optional) trained trace struct from a prior TRACE call;
%                    when provided, activates evaluation mode.
%
%   Outputs
%     out  - struct with fields:
%              space           convex-hull space footprint (measure, density, ...)
%              good{nalgos}    good-performance footprints
%              best{nalgos}    best-algorithm footprints
%              hard            beta-hard footprint (~beta instances)
%              summary         (nalgos+1 x 11) cell array performance table
%
%   TRACE3 reference
%     Simpson, D. et al. (2025). [TRACE3 paper citation.]
%
%   By: Mario Andres Munoz Acosta
%       School of Mathematics and Statistics
%       The University of Melbourne
%       Australia
%       2020 (TRACE3 refactor 2025)

narginchk(7, 8);
isEvalMode = nargin == 8;
is3D       = size(Z, 2) == 3;
nalgos     = size(Ybin, 2);
useLegacy  = isfield(opts, 'method') && strcmpi(opts.method, 'legacy');
if useLegacy && is3D
    warning('ISA:TRACE:legacyNo3D', ...
        'Legacy TRACE does not support 3D instance spaces; switching to TRACE3.');
    useLegacy = false;
end
% opts.pythiaSkip is set by buildIS.m when opts.pythia.skip=true. It must be
% checked explicitly: PYTHIA's skip-mode Yhat is a full-shape false(...)
% placeholder (not []), so isempty(Yhat) alone never detects the skip case.
pythiaSkipped   = isfield(opts, 'pythiaSkip') && opts.pythiaSkip;
pythiaAvailable = ~isempty(Yhat) && ~pythiaSkipped;

% -------------------------------------------------------------------------
% Space measure from convex hull of all instances (computed once).
% Passed to TRACEbuild3 so stopping criteria are relative to bounded space.
if is3D
    [~, spaceArea] = convhull(Z);
else
    [~, spaceArea] = convhull(Z(:,1), Z(:,2));
end
if is3D; measureLabel = 'Volume'; else; measureLabel = 'Area'; end

% =========================================================================
% EVALUATION MODE: re-score trained footprints on new instances
% =========================================================================
if isEvalMode
    fprintf('  -> TRACE is evaluating trained footprints on new instances.\n');
    out = trainedTrace;
    % Normalise backward-compat: old models stored .area instead of .measure
    if ~isfield(trainedTrace.space, 'measure')
        trainedTrace.space.measure  = trainedTrace.space.area;
        trainedTrace.space.density  = trainedTrace.space.density;
    end
    if ~isfield(trainedTrace.space, 'measureLabel')
        trainedTrace.space.measureLabel = measureLabel;
    end
    out.space = trainedTrace.space;
    ngood = numel(trainedTrace.good);
    nbest = numel(trainedTrace.best);
    for i = 1:min(nalgos, ngood)
        out.good{i} = TRACErescore(trainedTrace.good{i}, Z, Ybin(:,i), is3D);
    end
    for i = ngood+1:nalgos
        out.good{i} = TRACEthrow3(is3D);
    end
    for i = 1:min(nalgos, nbest)
        out.best{i} = TRACErescore(trainedTrace.best{i}, Z, P==i, is3D);
    end
    for i = nbest+1:nalgos
        out.best{i} = TRACEthrow3(is3D);
    end
    out.hard    = TRACErescore(trainedTrace.hard, Z, ~beta, is3D);
    out.summary = TRACEsummaryTable(out.good, out.best, algolabels, trainedTrace.space);
    fprintf('  -> Evaluation complete.\n');
    return;
end

% =========================================================================
% TRAINING MODE — LEGACY
% =========================================================================
if useLegacy
    fprintf('-------------------------------------------------------------------------\n');
    fprintf('  -> TRACE (legacy) is building footprints.\n');
    % Contradiction removal defaults to true for legacy (spec 4.1)
    useContra = ~isfield(opts, 'contra') || opts.contra;
    out = TRACE_legacy(Z, Ybin, P, beta, algolabels, opts, useContra);
    out = normalizeLegacyOut(out, measureLabel, nalgos);
    out.summary = TRACEsummaryTable(out.good, out.best, algolabels, out.space);
    fprintf('  -> TRACE (legacy) has completed. Footprint analysis results:\n\n');
    disp(out.summary);
    return;
end

% =========================================================================
% TRAINING MODE — TRACE3 (default)
% =========================================================================
fprintf('-------------------------------------------------------------------------\n');
fprintf('  -> TRACE3 is calculating the space %s and density.\n', lower(measureLabel));
out.space.measure      = spaceArea;
out.space.measureLabel = measureLabel;
out.space.elements     = size(Z, 1);
out.space.density      = out.space.elements / spaceArea;
out.space.purity       = 1;
out.space.polygon      = [];
fprintf('    -> Space %s: %s | Space density: %s\n', ...
    measureLabel, num2str(spaceArea), num2str(out.space.density));

if ~pythiaAvailable
    warning('ISA:TRACE3:noPYTHIA', ...
        'PYTHIA predictions unavailable; using true labels only (Zu = {yi=1}).');
end

if exist('gcp', 'file') == 2
    pool = gcp('nocreate');
    nworkers = 0;
    if ~isempty(pool), nworkers = pool.NumWorkers; end
else
    nworkers = 0;
end

fprintf('-------------------------------------------------------------------------\n');
fprintf('  -> TRACE3 is calculating the algorithm footprints.\n');
good = cell(1, nalgos);
best = cell(1, nalgos);
parfor (i = 1:nalgos, nworkers)
    t = tic;
    yhat_i = [];
    if pythiaAvailable, yhat_i = Yhat(:,i); end
    fprintf('    -> Good performance footprint for ''%s''\n', algolabels{i});
    good{i} = TRACEbuild3(Z, Ybin(:,i), yhat_i, spaceArea, opts);
    fprintf('    -> Best performance footprint for ''%s''\n', algolabels{i});
    best{i} = TRACEbuild3(Z, P==i, yhat_i, spaceArea, opts);
    fprintf('    -> Algorithm ''%s'' completed. Elapsed time: %.2fs\n', algolabels{i}, toc(t));
end
out.good = good;
out.best = best;

fprintf('-------------------------------------------------------------------------\n');
fprintf('  -> TRACE3 is calculating the beta-footprint.\n');
out.hard = TRACEbuild3(Z, ~beta, [], spaceArea, opts);

fprintf('-------------------------------------------------------------------------\n');
fprintf('  -> Preparing the summary table.\n');
out.summary = TRACEsummaryTable(out.good, out.best, algolabels, out.space);
fprintf('  -> TRACE3 has completed. Footprint analysis results:\n\n');
disp(out.summary);

end

% =========================================================================
% SUBFUNCTIONS
% =========================================================================

function footprint = TRACEbuild3(Z, Ybin, Yhat, spaceArea, opts)
% Build a TRACE3 footprint for one binary performance vector.
% Yhat = [] triggers fallback: Zu = {zi | yi=1} with no prediction filter.
is3D = size(Z, 2) == 3;

% Step 3: Zu = {zi | yhat_i=1 AND ybin_i=1}
if isempty(Yhat)
    Zu = Z(logical(Ybin), :);
else
    Zu = Z(logical(Yhat) & logical(Ybin), :);
end
Zu = unique(Zu, 'rows');

if size(Zu, 1) <= opts.minInstances
    footprint = TRACEthrow3(is3D);
    return;
end

% Step 4: build alpha-shape at default (minimum enclosing) alpha
as = alphaShape(Zu);

% Step 5: compute initial metrics
[footprint, valid] = TRACEmetrics3(as, Z, Ybin, is3D);
if ~valid || footprint.measure < opts.minAreaFrac * spaceArea
    footprint = TRACEthrow3(is3D);
    return;
end
if footprint.purity >= opts.PI
    return;
end

% Steps 6-7: iterate through alpha spectrum to tighten purity
AS = alphaSpectrum(as);
if numel(AS) < 2
    return;
end
alphaVec = linspace(as.Alpha, min(AS), 101);
alphaVec = alphaVec(2:end);  % first entry already evaluated above
for ii = 1:numel(alphaVec)
    as.Alpha = alphaVec(ii);
    if is3D
        as.RegionThreshold = volume(as) / 20;
    else
        as.RegionThreshold = area(as) / 20;
    end
    [footprint, valid] = TRACEmetrics3(as, Z, Ybin, is3D);
    if ~valid || footprint.measure < opts.minAreaFrac * spaceArea
        footprint = TRACEthrow3(is3D);
        return;
    end
    if footprint.purity >= opts.PI
        return;
    end
end
% Alpha spectrum exhausted — return the best footprint found.
end

% =========================================================================
function [footprint, valid] = TRACEmetrics3(as, Z, Ybin, is3D)
% Compute footprint metrics from an alphaShape object.
valid = true;
footprint.polygon = as;
if is3D
    m = volume(as);
    footprint.measureLabel = 'Volume';
else
    m = area(as);
    footprint.measureLabel = 'Area';
end
if m <= 0 || isinf(as.Alpha)
    valid = false;
    footprint.measure      = 0;
    footprint.elements     = 0;
    footprint.goodElements = 0;
    footprint.density      = 0;
    footprint.purity       = 0;
    return;
end
footprint.measure      = m;
footprint.elements     = sum(inShape(as, Z));
footprint.goodElements = sum(inShape(as, Z(logical(Ybin),:)));
if footprint.elements == 0
    valid = false;
    footprint.density = 0;
    footprint.purity  = 0;
    return;
end
footprint.density = footprint.elements / m;
footprint.purity  = footprint.goodElements / footprint.elements;
end

% =========================================================================
function footprint = TRACEthrow3(is3D)
fprintf('        -> There are not enough instances to calculate a footprint.\n');
footprint.polygon      = [];
footprint.measure      = 0;
footprint.measureLabel = 'Area';
if is3D, footprint.measureLabel = 'Volume'; end
footprint.elements     = 0;
footprint.goodElements = 0;
footprint.density      = 0;
footprint.purity       = 0;
end

% =========================================================================
function footprint = TRACErescore(trainedFp, Z, Ybin, is3D)
% Re-evaluate a trained footprint against new instances Z.
footprint = trainedFp;
% Backward-compat: old models stored .area; normalise to .measure
if isfield(footprint, 'area') && ~isfield(footprint, 'measure')
    footprint.measure = footprint.area;
    if is3D; footprint.measureLabel = 'Volume'; else; footprint.measureLabel = 'Area'; end
end
poly = trainedFp.polygon;
if isempty(poly)
    footprint.elements     = 0;
    footprint.goodElements = 0;
    footprint.density      = 0;
    footprint.purity       = 0;
    return;
end
if isa(poly, 'alphaShape')
    inside = inShape(poly, Z);
elseif isa(poly, 'polyshape')
    inside = isinterior(poly, Z);
else
    footprint.elements = 0; footprint.goodElements = 0;
    footprint.density  = 0; footprint.purity       = 0;
    return;
end
footprint.elements     = sum(inside);
footprint.goodElements = sum(inside & logical(Ybin));
m = footprint.measure;
if footprint.elements == 0 || m == 0
    footprint.density = 0;
    footprint.purity  = 0;
else
    footprint.density = footprint.elements / m;
    footprint.purity  = footprint.goodElements / footprint.elements;
end
end

% =========================================================================
function summary = TRACEsummaryTable(good, best, algolabels, space)
nGood = length(good);
nBest = length(best);
nLabels = numel(algolabels);
if nGood ~= nBest || nGood ~= nLabels
    error('ISA:TRACE:summaryMismatch', ...
        'TRACEsummaryTable: good=%d, best=%d, labels=%d must all match.', ...
        nGood, nBest, nLabels);
end
nalgos = nGood;
ml = space.measureLabel;
summary = cell(nalgos+1, 11);
summary(1, 2:end) = {[ml '_Good'],            [ml '_Good_Normalized'], ...
                     'Density_Good',           'Density_Good_Normalized', ...
                     'Purity_Good', ...
                     [ml '_Best'],             [ml '_Best_Normalized'], ...
                     'Density_Best',           'Density_Best_Normalized', ...
                     'Purity_Best'};
summary(2:end, 1) = algolabels;
for i = 1:nalgos
    row = [TRACEsummaryRow(good{i}, space.measure, space.density), ...
           TRACEsummaryRow(best{i}, space.measure, space.density)];
    summary(i+1, 2:end) = num2cell(round(row, 3));
end
end

% =========================================================================
function row = TRACEsummaryRow(fp, spaceMeasure, spaceDensity)
m = fp.measure;
row = [m, m/spaceMeasure, fp.density, fp.density/spaceDensity, fp.purity];
row(isnan(row)) = 0;
end

% =========================================================================
function out = normalizeLegacyOut(out, measureLabel, nalgos)
% Add .measure / .measureLabel to all footprints from a legacy out struct.
out.space = addMeasureField(out.space, measureLabel);
for i = 1:nalgos
    out.good{i} = addMeasureField(out.good{i}, measureLabel);
    out.best{i} = addMeasureField(out.best{i}, measureLabel);
end
out.hard = addMeasureField(out.hard, measureLabel);
end

% =========================================================================
function fp = addMeasureField(fp, ml)
if isfield(fp, 'area')
    fp.measure = fp.area;
elseif ~isfield(fp, 'measure')
    fp.measure = 0;
end
fp.measureLabel = ml;
end
