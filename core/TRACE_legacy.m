function out = TRACE_legacy(Z, Ybin, P, beta, algolabels, opts, useContra)
% TRACE_legacy  Legacy DBSCAN+polyshape footprint algorithm.
%
%   out = TRACE_legacy(Z, Ybin, P, beta, algolabels, opts, useContra)
%
%   Selectable via opts.trace.method = 'legacy' in TRACE.m. Uses DBSCAN
%   clustering followed by boundary/polyshape construction per cluster.
%   Contradiction removal between best-algorithm footprints is applied when
%   useContra = true (the default for legacy mode).
%
%   Inputs / outputs mirror the old TRACE.m interface. Footprint structs use
%   .area (not .measure); TRACE.m normalises naming after calling this function.

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

if nargin < 7, useContra = true; end

if exist('gcp', 'file') == 2
    pool = gcp('nocreate');
    nworkers = 0;
    if ~isempty(pool), nworkers = pool.NumWorkers; end
else
    nworkers = 0;
end

ninst  = size(Z, 1);
nalgos = size(Ybin, 2);

% Space footprint (all instances)
fprintf('[TRACE] TRACE (legacy) is calculating the space area and density.\n');
out.space = TRACEbuild(Z, true(ninst,1), opts);
fprintf('[TRACE] Space area: %s | Space density: %s\n', ...
    num2str(out.space.area), num2str(out.space.density));

% Per-algorithm footprints
fprintf('[TRACE] TRACE (legacy) is calculating the algorithm footprints.\n');
good = cell(1, nalgos);
best = cell(1, nalgos);
parfor (i = 1:nalgos, nworkers)
    t = tic;
    fprintf('[TRACE] Good performance footprint for ''%s''\n', algolabels{i});
    good{i} = TRACEbuild(Z, Ybin(:,i), opts);
    fprintf('[TRACE] Best performance footprint for ''%s''\n', algolabels{i});
    best{i} = TRACEbuild(Z, P==i, opts);
    fprintf('[TRACE] Algorithm ''%s'' completed. Elapsed time: %.2fs\n', algolabels{i}, toc(t));
end
out.good = good;
out.best = best;

% Contradiction removal
if useContra
    fprintf('[TRACE] TRACE (legacy) is detecting and removing contradictory footprint sections.\n');
    for i = 1:nalgos
        fprintf('[TRACE] Base algorithm ''%s''\n', algolabels{i});
        tBase = tic;
        for j = i+1:nalgos
            fprintf('[TRACE] Comparing ''%s'' with ''%s''\n', algolabels{i}, algolabels{j});
            tTest = tic;
            [out.best{i}, out.best{j}] = TRACEcontra(out.best{i}, out.best{j}, ...
                                                      Z, P==i, P==j, opts);
            fprintf('[TRACE] Completed. Elapsed time: %.2fs\n', toc(tTest));
        end
        fprintf('[TRACE] Base algorithm ''%s'' completed. Elapsed time: %.2fs\n', algolabels{i}, toc(tBase));
    end
end

% Beta hard footprint
fprintf('[TRACE] TRACE (legacy) is calculating the beta-footprint.\n');
out.hard = TRACEbuild(Z, ~beta, opts);

end

% =========================================================================
% SUBFUNCTIONS
% =========================================================================

function footprint = TRACEbuild(Z, Ybin, opts)
Ig = unique(Z(logical(Ybin),:), 'rows');
if size(Ig, 1) < 3
    footprint = TRACEthrow;
    return;
else
    footprint = struct;
end

nn    = max(min(ceil(sum(Ybin)/20), 50), 3);
class = TRACEdbscan(Ig, nn);
flag  = false;
for i = 1:max(class)
    polydata = Ig(class==i, :);
    polydata = polydata(boundary(polydata,1), :);
    aux = TRACEfitpoly(polydata, Z, Ybin, opts);
    if ~isempty(aux)
        if ~flag
            footprint.polygon = aux;
            flag = true;
        else
            footprint.polygon = union(footprint.polygon, aux);
        end
    end
end
if isfield(footprint, 'polygon') && ~isempty(footprint.polygon)
    footprint.polygon      = rmslivers(footprint.polygon, 1e-2);
    footprint.area         = area(footprint.polygon);
    footprint.elements     = sum(isinterior(footprint.polygon, Z));
    footprint.goodElements = sum(isinterior(footprint.polygon, Z(logical(Ybin),:)));
    footprint.density      = footprint.elements / footprint.area;
    footprint.purity       = footprint.goodElements / footprint.elements;
else
    footprint = TRACEthrow;
end
end

% =========================================================================
function [base, test] = TRACEcontra(base, test, Z, Ybase, Ytest, opts)
if isempty(base.polygon) || isempty(test.polygon)
    return;
end

maxtries = 3;
numtries = 1;
contradiction = intersect(base.polygon, test.polygon);
while contradiction.NumRegions ~= 0 && numtries <= maxtries
    numElements          = sum(isinterior(contradiction, Z));
    numGoodElementsBase  = sum(isinterior(contradiction, Z(logical(Ybase),:)));
    numGoodElementsTest  = sum(isinterior(contradiction, Z(logical(Ytest),:)));
    purityBase = numGoodElementsBase / numElements;
    purityTest = numGoodElementsTest / numElements;
    if purityBase > purityTest
        carea = area(contradiction) / area(test.polygon);
        fprintf('[TRACE] %.1f%% of the test footprint is contradictory.\n', round(100.*carea,1));
        test.polygon = subtract(test.polygon, contradiction);
        if numtries < maxtries
            test.polygon = TRACEtight(test.polygon, Z, Ytest, opts);
        end
    elseif purityTest > purityBase
        carea = area(contradiction) / area(base.polygon);
        fprintf('[TRACE] %.1f%% of the base footprint is contradictory.\n', round(100.*carea,1));
        base.polygon = subtract(base.polygon, contradiction);
        if numtries < maxtries
            base.polygon = TRACEtight(base.polygon, Z, Ybase, opts);
        end
    else
        fprintf('[TRACE] Purity of the contradicting areas is equal for both footprints.\n');
        fprintf('[TRACE] Ignoring the contradicting area.\n');
        break;
    end
    if isempty(base.polygon) || isempty(test.polygon)
        break;
    else
        contradiction = intersect(base.polygon, test.polygon);
    end
    numtries = numtries + 1;
end

if isempty(base.polygon)
    base = TRACEthrow;
else
    base.area         = area(base.polygon);
    base.elements     = sum(isinterior(base.polygon, Z));
    base.goodElements = sum(isinterior(base.polygon, Z(logical(Ybase),:)));
    base.density      = base.elements / base.area;
    base.purity       = base.goodElements / base.elements;
end
if isempty(test.polygon)
    test = TRACEthrow;
else
    test.area         = area(test.polygon);
    test.elements     = sum(isinterior(test.polygon, Z));
    test.goodElements = sum(isinterior(test.polygon, Z(logical(Ytest),:)));
    test.density      = test.elements / test.area;
    test.purity       = test.goodElements / test.elements;
end
end

% =========================================================================
function polygon = TRACEtight(polygon, Z, Ybin, opts)
splits   = regions(polygon);
nregions = length(splits);
flags    = true(1, nregions);
for i = 1:nregions
    criteria = isinterior(splits(i), Z) & logical(Ybin);
    polydata = Z(criteria, :);
    if size(polydata,1) < 3
        flags(i) = false;
        continue;
    end
    aux = TRACEfitpoly(polydata(boundary(polydata,1),:), Z, Ybin, opts);
    if isempty(aux)
        flags(i) = false;
        continue;
    end
    splits(i) = aux;
end
if any(flags)
    polygon = union(splits(flags));
else
    polygon = [];
end
end

% =========================================================================
function polygon = TRACEfitpoly(polydata, Z, Ybin, opts)
warning('off', 'MATLAB:polyshape:repairedBySimplify');
if size(polydata,1) < 3
    polygon = [];
    warning('on', 'MATLAB:polyshape:repairedBySimplify');
    return;
end
polygon = polyshape(polydata, 'Simplify', true);
polygon = rmslivers(polygon, 5e-2);
if ~all(Ybin)
    if polygon.NumRegions < 1
        polygon = [];
        warning('on', 'MATLAB:polyshape:repairedBySimplify');
        return;
    end
    tri  = triangulation(polygon);
    nrow = size(tri.ConnectivityList, 1);
    for ii = 1:nrow
        tridata  = tri.Points(tri.ConnectivityList(ii,:), :);
        piece    = polyshape(tridata, 'Simplify', true);
        elements = sum(isinterior(piece, Z));
        goodElements = sum(isinterior(piece, Z(logical(Ybin),:)));
        if elements == 0 || opts.PI > (goodElements/elements)
            polygon = subtract(polygon, piece);
        end
    end
end
warning('on', 'MATLAB:polyshape:repairedBySimplify');
end

% =========================================================================
function footprint = TRACEthrow
fprintf('[TRACE] There are not enough instances to calculate a footprint.\n');
footprint.polygon      = [];
footprint.area         = 0;
footprint.elements     = 0;
footprint.goodElements = 0;
footprint.density      = 0;
footprint.purity       = 0;
end

% =========================================================================
% DBSCAN implementation (Daszykowski 2004). Renamed from dbscan() to
% TRACEdbscan() to avoid shadowing the toolbox function added in R2019a.
% =========================================================================
function [class, type] = TRACEdbscan(x, k, Eps)
m = size(x,1);
if nargin < 3 || isempty(Eps)
    Eps = TRACEepsilon(x, k);
end
x       = [(1:m)' x];
[m, n]  = size(x);
type    = zeros(1, m);
no      = 1;
touched = zeros(m, 1);
class   = zeros(1, m);
for i = 1:m
    if touched(i) == 0
        ob  = x(i,:);
        D   = TRACEdist(ob(2:n), x(:,2:n));
        ind = find(D <= Eps);
        if length(ind) > 1 && length(ind) < k+1
            type(i) = 0; class(i) = 0;
        end
        if length(ind) == 1
            type(i) = -1; class(i) = -1; touched(i) = 1;
        end
        if length(ind) >= k+1
            type(i)     = 1;
            class(ind)  = ones(length(ind),1) * max(no);
            while ~isempty(ind)
                ob  = x(ind(1),:);
                touched(ind(1)) = 1;
                ind(1) = [];
                D  = TRACEdist(ob(2:n), x(:,2:n));
                i1 = find(D <= Eps);
                if length(i1) > 1
                    class(i1) = no;
                    if length(i1) >= k+1
                        type(ob(1)) = 1;
                    else
                        type(ob(1)) = 0;
                    end
                    for j = 1:length(i1)
                        if touched(i1(j)) == 0
                            touched(i1(j)) = 1;
                            ind = [ind i1(j)]; %#ok<AGROW>
                            class(i1(j)) = no;
                        end
                    end
                end
            end
            no = no + 1;
        end
    end
end
i1 = find(class == 0);
class(i1) = -1;
type(i1)  = -1;
end

% =========================================================================
function Eps = TRACEepsilon(x, k)
[m, n] = size(x);
Eps = ((prod(max(x)-min(x)) * k * gamma(.5*n+1)) / (m * sqrt(pi.^n))) .^ (1/n);
end

% =========================================================================
function D = TRACEdist(i, x)
[m, n] = size(x);
D = sqrt(sum((((ones(m,1)*i) - x).^2)'));
if n == 1
    D = abs((ones(m,1)*i - x))';
end
end
