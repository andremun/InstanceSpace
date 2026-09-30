classdef OptionsReviewTest < matlab.unittest.TestCase
% OptionsReviewTest  Regression cases for the September 2026 review.
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


    methods (Test)
        function testGroupJsonRoundTrip(tc)
            groups = {{[1 2],[3 4]}, {[1 2],3}, {1,2}, {}};
            for i = 1:numel(groups)
                opts.pilot.viewGroups = groups{i};
                decoded = ISAvalidateOpts(jsondecode(jsonencode(opts)));
                tc.verifyEqual(cellfun(@(g) g(:)',decoded.pilot.viewGroups(:),'UniformOutput',false), ...
                    cellfun(@(g) g(:)',groups{i}(:),'UniformOutput',false));
            end
        end
        function testCameraDirection(tc)
            scriptfcn;
            fig = figure('Visible','off'); tc.addTeardown(@() close(fig));
            for direction = [1 0 0;0 1 0;1 2 3]'
                [az,el] = cart2sph(direction(1),direction(2),direction(3));
                v = struct('azimuth',az,'elevation',el,'groups',{{1}});
                angles = resolveViewAngle(v,1);
                drawScatter([0 0 0;1 2 3],[0;1],'test',angles);
                ax = gca; actual = ax.CameraPosition-ax.CameraTarget;
                tc.verifyEqual(actual/norm(actual),direction'/norm(direction),'AbsTol',1e-5);
            end
        end
    end
end
