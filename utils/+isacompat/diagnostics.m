function report = diagnostics()
% diagnostics  Inspect the runtime and API locations without changing paths.
% Installed Octave packages are recorded but never loaded or installed.
% API presence is not a compatibility guarantee. Run portable smoke tests
% separately; see docs/octave-compatibility.md for the serial release scope.
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

report.schemaVersion = 1;
report.isOctave = isacompat.isOctave();
report.version = version;
report.engine = 'MATLAB';
report.packages = struct('name', {}, 'version', {}, 'loaded', {});
if report.isOctave
    report.engine = 'Octave';
    installed = pkg('list');
    for k = 1:numel(installed)
        p = installed{k};
        report.packages(k) = struct('name', p.name, ...
            'version', p.version, 'loaded', logical(p.loaded));
    end
end
names = {'readtable','table','string','nanmedian','boxcox','zscore', ...
    'corr','pdist','plsregress','optimoptions','fminunc','ga','kmeans', ...
    'cvpartition','fitcknn','fitcsvm','fitctree','fitcnb','fitclinear', ...
    'fitcensemble','fitSVMPosterior','templateTree','sobolset','scramble','bayesopt', ...
    'alphaShape','polyshape','exportgraphics','gcp','rng'};
report.apis = struct();
for k = 1:numel(names)
    location = which(names{k});
    report.apis.(names{k}) = struct('available', ~isempty(location), ...
                                  'location', location);
end
report.octaveScope = 'Serial 2D/3D build/save/load/explore; KNN none/Sobol; TRACE3; CSV/PNG (Qt)';
end
