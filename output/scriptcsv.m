function scriptcsv(container,rootdir)
% scriptcsv  Write a model or explore() result to CSV files in rootdir.
%
%   scriptcsv(container,rootdir)
%
%   Writes projected coordinates, feature/performance tables, algorithm
%   selections, and footprint boundary points, sized for 2D or 3D
%   projections according to size(container.pilot.Z,2).
%
%   Inputs
%     container - model struct from buildIS/InstanceSpace.build() or a testResults entry from exploreIS/InstanceSpace.explore().
%     rootdir - destination directory (trailing slash required).
%   Outputs
%     none - writes CSV files to rootdir as a side effect.

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

scriptfcn;

nalgos = size(container.data.Y,2);
fprintf('[OUTPUT] Writing the data on CSV files for posterior analysis.\n');
% -------------------------------------------------------------------------
% Determine dimensionality
ndim = size(container.pilot.Z, 2);
if ndim == 3
    zcols = {'z_1','z_2','z_3'};
else
    zcols = {'z_1','z_2'};
end

% These names belong to the toolkit's geometry export namespace. Remove
% prior geometry, including algorithms absent from the current portfolio.
old = dir(fullfile(rootdir, '*.csv'));
for j = 1:numel(old)
    name = old(j).name;
    if ~isempty(regexp(name, '^footprint_.+_(best|good)(_faces)?\.csv$', 'once')) || ...
            ismember(name, {'bounds.csv','bounds_prunned.csv','bounds_faces.csv','bounds_prunned_faces.csv'})
        delete(fullfile(rootdir,name));
    end
end
geometry = struct('algorithm',{},'kind',{},'vertices',{},'faces',{},'status',{});
for i = 1:nalgos
    for kind = {'best','good'}
        fp = container.trace.(kind{1}){i};
        faces = zeros(0,3);
        if ndim == 3 && isfield(fp,'polygon') && isacompat.isAlphaShape(fp.polygon) && ~isempty(fp.polygon.Points)
            [faces,verts] = boundaryFacets(fp.polygon);
        else
            verts = footprintBoundary(fp);
        end
        if isempty(verts), verts = zeros(0,ndim); end
        base = ['footprint_' container.data.algolabels{i} '_' kind{1}];
        writeArray2CSV(verts,zcols,makeBndLabels(verts),fullfile(rootdir,[base '.csv']));
        facefile = '';
        if ndim == 3
            facefile = [base '_faces.csv'];
            writematrixWithHeaders(faces,fullfile(rootdir,facefile));
        end
        status = 'available';
        if isempty(verts), status = 'empty'; end
        geometry(end+1) = struct('algorithm',container.data.algolabels{i}, ...
            'kind',kind{1},'vertices',[base '.csv'],'faces',facefile,'status',status); %#ok<AGROW>
    end
end
manifest = struct('dimensions',ndim,'indexBase',1, ...
    'rings','2D boundary rings are NaN-separated. Nested rings represent holes.', ...
    'footprints',geometry);
fid = fopen(fullfile(rootdir,'geometry_manifest.json'),'w');
if fid < 0, error('ISA:scriptcsv:manifestWrite','Cannot write geometry manifest.'); end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'%s',jsonencode(manifest));
clear cleanup;

writeArray2CSV(container.pilot.Z, zcols, ...
               container.data.instlabels, ...
               [rootdir 'coordinates.csv']);
if isfield(container,'cloist')
    if ndim == 3
        if ~isfield(container.cloist,'ZedgeFaces') || ~isfield(container.cloist,'ZecorrFaces')
            error('ISA:scriptcsv:missingFaces','Rebuild CLOISTER before exporting this legacy 3D boundary.');
        end
        writematrixWithHeaders(container.cloist.ZedgeFaces,fullfile(rootdir,'bounds_faces.csv'));
        writematrixWithHeaders(container.cloist.ZecorrFaces,fullfile(rootdir,'bounds_prunned_faces.csv'));
    end
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
              [rootdir 'classifier_table.csv']);
end

function writematrixWithHeaders(faces,file)
% Face indices refer to the corresponding vertex CSV rows, starting at one.
writetable(array2table(faces,'VariableNames',{'vertex_1','vertex_2','vertex_3'}),file);
end
