classdef OctaveCompatibilityTest < matlab.unittest.TestCase
% OctaveCompatibilityTest  Run the same portable contracts in MATLAB CI.
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
        function sharedGeometryContracts(~)
            folder = fullfile(fileparts(mfilename('fullpath')), 'portable');
            oldPath = path;
            cleanup = onCleanup(@() path(oldPath)); %#ok<NASGU>
            addpath(folder);
            runGeometrySmoke();
        end

        function literalArchiveStructs(~)
            folder = fullfile(fileparts(mfilename('fullpath')), 'portable');
            oldPath = path;
            cleanup = onCleanup(@() path(oldPath)); %#ok<NASGU>
            addpath(folder);
            runArchiveSmoke();
        end

        function sharedOptimizationContracts(~)
            folder = fullfile(fileparts(mfilename('fullpath')), 'portable');
            oldPath = path;
            cleanup = onCleanup(@() path(oldPath)); %#ok<NASGU>
            addpath(folder);
            runOptimizationSmoke();
        end

        function sharedLearningContracts(~)
            folder = fullfile(fileparts(mfilename('fullpath')), 'portable');
            oldPath = path;
            cleanup = onCleanup(@() path(oldPath)); %#ok<NASGU>
            addpath(folder);
            runLearningSmoke();
        end
        function sharedDeterministicContracts(tc)
            folder = fullfile(fileparts(mfilename('fullpath')), 'portable');
            oldPath = path;
            cleanup = onCleanup(@() path(oldPath)); %#ok<NASGU>
            addpath(folder);
            report = runPortableSmoke();
            tc.verifyFalse(report.isOctave);
        end
    end
end
