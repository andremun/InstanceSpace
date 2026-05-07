function scriptcsv(container,rootdir)
% -------------------------------------------------------------------------
% csvscript.m
% -------------------------------------------------------------------------
%
% By: Mario Andres Munoz Acosta
%     School of Mathematics and Statistics
%     The University of Melbourne
%     Australia
%     2020
%
% -------------------------------------------------------------------------

scriptfcn;

nalgos = size(container.data.Y,2);
ndim = size(container.pilot.Z, 2);
if ndim == 3
    zcols = {'z_1','z_2','z_3'};
else
    zcols = {'z_1','z_2'};
end
fprintf('=========================================================================\n');
fprintf('-> Writing the data on CSV files for posterior analysis.\n');
% -------------------------------------------------------------------------
for i=1:nalgos
    verts = footprintBoundary(container.trace.best{i});
    if ~isempty(verts)
        writeArray2CSV(verts, zcols, makeBndLabels(verts), ...
                       [rootdir 'footprint_' container.data.algolabels{i} '_best.csv']);
    end
    verts = footprintBoundary(container.trace.good{i});
    if ~isempty(verts)
        writeArray2CSV(verts, zcols, makeBndLabels(verts), ...
                       [rootdir 'footprint_' container.data.algolabels{i} '_good.csv']);
    end
%     if isfield(container.trace.bad{i},'polygon') && ~isempty(container.trace.bad{i}.polygon)
%         writeArray2CSV(container.trace.bad{i}.polygon.Vertices, ...
%                        {'z_1','z_2'},...
%                        makeBndLabels(container.trace.bad{i}.polygon.Vertices),...
%                        [rootdir 'footprint_' container.data.algolabels{i} '_bad.csv']);
%     end
end

writeArray2CSV(container.pilot.Z, zcols, ...
               container.data.instlabels, ...
               [rootdir 'coordinates.csv']);
if isfield(container,'cloist')
    writeArray2CSV(container.cloist.Zedge, zcols, ...
                   makeBndLabels(container.cloist.Zedge), ...
                   [rootdir 'bounds.csv']);
    writeArray2CSV(container.cloist.Zecorr, zcols, ...
                   makeBndLabels(container.cloist.Zecorr), ...
                   [rootdir 'bounds_prunned.csv']);
end
writeArray2CSV(container.data.Xraw(:, container.featsel.idx), ...
               container.data.featlabels, ...
               container.data.instlabels, ...
               [rootdir 'feature_raw.csv']);
writeArray2CSV(container.data.X, ...
               container.data.featlabels, ...
               container.data.instlabels, ...
               [rootdir 'feature_process.csv']);
writeArray2CSV(container.data.Yraw, ...
               container.data.algolabels, ...
               container.data.instlabels, ...
               [rootdir 'algorithm_raw.csv']);
writeArray2CSV(container.data.Y, ...
               container.data.algolabels, ...
               container.data.instlabels, ...
               [rootdir 'algorithm_process.csv']);
writeArray2CSV(container.data.Ybin, ...
               container.data.algolabels, ...
               container.data.instlabels, ...
               [rootdir 'algorithm_bin.csv']);
writeArray2CSV(container.data.numGoodAlgos, {'NumGoodAlgos'}, ...
               container.data.instlabels, ...
               [rootdir 'good_algos.csv']);
writeArray2CSV(container.data.beta, {'IsBetaEasy'}, ...
               container.data.instlabels, ...
               [rootdir 'beta_easy.csv']);
writeArray2CSV(container.data.P, {'Best_Algorithm'}, ...
               container.data.instlabels, ...
               [rootdir 'portfolio.csv']);
writeArray2CSV(container.pythia.Yhat, ...
               container.data.algolabels, ...
               container.data.instlabels, ...
               [rootdir 'algorithm_svm.csv']);
writeArray2CSV(container.pythia.selection0, ...
               {'Best_Algorithm'}, ...
               container.data.instlabels, ...
               [rootdir 'portfolio_svm.csv']);
writeCell2CSV(container.trace.summary(2:end,[3 5 6 8 10 11]), ...
              container.trace.summary(1,[3 5 6 8 10 11]),...
              container.trace.summary(2:end,1),...
              [rootdir 'footprint_performance.csv']);
if isfield(container.pilot,'summary')
    writeCell2CSV(container.pilot.summary(2:end,2:end), ...
                  container.pilot.summary(1,2:end),...
                  container.pilot.summary(2:end,1), ...
                  [rootdir 'projection_matrix.csv']);
end
writeCell2CSV(container.pythia.summary(2:end,2:end), ...
              container.pythia.summary(1,2:end), ...
              container.pythia.summary(2:end,1), ...
              [rootdir 'svm_table.csv']);
end

% =========================================================================
function verts = footprintBoundary(fp)
% Extract ordered 2-D boundary vertices from a footprint struct.
% Handles both polyshape (legacy) and alphaShape (TRACE3) polygons.
verts = [];
if ~isfield(fp, 'polygon') || isempty(fp.polygon)
    return;
end
poly = fp.polygon;
if isa(poly, 'alphaShape')
    if size(poly.Points, 2) == 3
        % 3D boundary CSV export not supported; skip silently.
        % Full 3D boundary extraction is deferred to a later phase.
        return;
    end
    [bf, bv] = boundaryFacets(poly);
    if isempty(bf)
        return;
    end
    verts = traceAlphaBoundary(bf, bv);
elseif isa(poly, 'polyshape')
    verts = poly.Vertices;
end
end

% =========================================================================
function verts = traceAlphaBoundary(bf, bv)
% Trace an ordered closed polygon from a 2-D boundary-facets edge list.
% bf: (m x 2) edge index pairs into bv rows; bv: (n x 2) point coordinates.
% Works correctly for simple (single-region) connected alpha shapes.
n = size(bv, 1);
if n == 0, verts = []; return; end
% Build adjacency: each vertex has exactly 2 neighbours on the boundary.
adj = zeros(n, 2);
cnt = zeros(n, 1);
for k = 1:size(bf, 1)
    v1 = bf(k,1); v2 = bf(k,2);
    cnt(v1) = cnt(v1)+1;
    if cnt(v1) <= 2, adj(v1, cnt(v1)) = v2; end
    cnt(v2) = cnt(v2)+1;
    if cnt(v2) <= 2, adj(v2, cnt(v2)) = v1; end
end
% Trace from vertex 1 following the boundary.
order = zeros(n, 1);
order(1) = 1;
prev = 0; curr = 1;
for k = 2:n
    nxt = adj(curr, adj(curr,:) ~= prev & adj(curr,:) ~= 0);
    if isempty(nxt), break; end
    order(k) = nxt(1);
    prev = curr;
    curr = order(k);
end
valid = order ~= 0;
verts = bv(order(valid), :);
end
