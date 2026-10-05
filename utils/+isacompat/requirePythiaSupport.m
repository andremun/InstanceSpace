function requirePythiaSupport(opts, trained)
% requirePythiaSupport  Bound the experimental Octave classifier surface.
% Called after PYTHIA fills its defaults. Does not load/install dependencies.
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

if ~isacompat.isOctave(), return; end
if compare_versions(version,'11.1.0','<')
    error('ISA:compat:runtimeVersion','Experimental Octave PYTHIA requires Octave >=11.1.0 and Statistics 2.0.0.');
end
required = {'fitcknn','cvpartition','nanmean','nanstd'};
for k = 1:numel(required)
    if isempty(which(required{k}))
        error('ISA:compat:missingPackages', ...
            'Missing %s. Load the tested Statistics 2.0.0 / Datatypes 1.5.0 environment before PYTHIA.',required{k});
    end
end
if nargin == 2
    if ~isfield(trained,'classifierType') || ~any(strcmpi(trained.classifierType,{'knn','none'}))
        error('ISA:compat:unsupportedFeature', ...
            'Octave evaluation currently supports KNN or skip-mode results only.');
    end
elseif ~opts.skip
    if ~strcmpi(opts.classifier,'knn') || ~any(strcmpi(opts.tuning,{'none','sobol'}))
        error('ISA:compat:unsupportedFeature', ...
            ['Experimental Octave PYTHIA supports classifier=''knn'', tuning=''none'' or ''sobol'' ' ...
             '(explicit opts.params for none). Other classifiers/search modes await validation.']);
    end
    if isfield(opts,'parallel') && opts.parallel
        isacompat.requireFeature('parallel');
    end
end
end
