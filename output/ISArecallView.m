function ISArecallView(fig, groupIdx)
% ISArecallView  Snap an open 3D instance-space figure to its stored
% PILOTviewpoint camera angle for a given algorithm group.
%
%   ISArecallView(fig, groupIdx)
%   ISArecallView(fig)          % use the default/global viewpoint
%
%   fig      - handle to a figure produced by scriptpng.m, or a .fig file
%              it wrote reopened via openfig/uiopen -- the viewpoint is
%              stored in the figure's UserData (scriptpng.m sets
%              fig.UserData.isaViewpoint), so it survives a save/load
%              round-trip through a .fig file even in a different MATLAB
%              session.
%   groupIdx - algorithm column index (as used in opts.pilot.viewGroups)
%              whose stored viewpoint to apply. Omit or pass [] for the
%              default/global viewpoint (the one feature/portfolio-level
%              plots use).
%
%   Useful after manually rotating a 3D footprint .fig while exploring it
%   interactively, to return to the optimised viewpoint (spec §5.2)
%   without having to recompute or look it up by hand.
%
%   Example:
%     fig = openfig('footprint_KNN.fig');
%     ISArecallView(fig, 3);   % back to the stored view for algorithm 3

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

if nargin < 2
    groupIdx = [];
end
if ~(isscalar(fig) && isgraphics(fig, 'figure'))
    error('ISA:ISArecallView:notAFigure', 'fig must be a valid figure handle.');
end
if ~isstruct(fig.UserData) || ~isfield(fig.UserData, 'isaViewpoint') || ...
        isempty(fig.UserData.isaViewpoint)
    error('ISA:ISArecallView:noStoredViewpoint', ...
        ['This figure has no stored PILOTviewpoint data (fig.UserData.isaViewpoint). ' ...
         'Only 3D figures written by scriptpng.m with a computed viewpoint carry this -- ' ...
         '2D projections, figures with no viewpoint, and figures from other sources do not.']);
end
viewpoint = fig.UserData.isaViewpoint;

% resolveViewAngle is scriptfcn.m's shared helper (also used by scriptpng.m
% itself); scriptfcn assigns it into ITS caller's workspace via assignin,
% so calling it here makes it available in this function's own workspace.
scriptfcn;
viewAngle = resolveViewAngle(viewpoint, groupIdx);
if isempty(viewAngle)
    error('ISA:ISArecallView:noViewAngle', ...
        'Could not resolve a stored viewpoint for the requested algorithm group.');
end

% findall(...,'Type','axes') also returns the internal axes MATLAB uses
% to implement legends/colorbars (Tag 'legend'/'Colorbar'), not just the
% main plot axes; applying view() to one of those instead would leave the
% actual footprint plot unrotated. Exclude them and apply the view to
% every remaining (real, data) axes rather than assuming there's only one.
ax = findall(fig, 'Type', 'axes');
tags = arrayfun(@(a) string(get(a, 'Tag')), ax);
ax = ax(~ismember(tags, ["legend", "Colorbar"]));
if isempty(ax)
    error('ISA:ISArecallView:noAxes', 'fig has no plot axes to apply the view to.');
end
for i = 1:numel(ax)
    view(ax(i), viewAngle(1), viewAngle(2));
end
end
