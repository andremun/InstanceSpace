% test_integration.m - ISA Toolkit full-pipeline regression suite: runs the
% full pipeline on the reference dataset and checks for numerical regressions.
%
% This is the exhaustive option-coverage suite, not a getting-started
% example -- for that, see example.m, which runs a single bare-bones
% configuration with a short guide to the handful of settings most users
% adjust first. This file exists to catch regressions across the option
% surface: new features and bug fixes should get a test method under
% tests/, not a new block appended to this file.
%
% #39: this used to be a bare script with hand-rolled try/catch bookkeeping
% per case. It's now a thin runner over matlab.unittest.TestCase classes in
% tests/ (PipelineOptionsTest, ClassApiTest, MigrationTest, RegressionTest),
% using MATLAB's own test framework for assertions, parameterisation
% (PipelineOptionsTest.OptionCase), and reporting -- run via runtests()
% would work too, but this file additionally wires in the code coverage
% plugin and preserves the EOF:SUCCESS/EOF:ERROR sentinel convention
% external automation may still scrape for (buildIS.m/exploreIS.m use the
% same convention independently of this file, so it isn't purely internal
% to what's being migrated here).
%
% IMPORTANT: options.json is NOT the source of truth for these settings.
% Each test class writes a fresh options.json into its own case
% subdirectory (test/data/<case_name>/options.json) on every run, from the
% opts struct built in MATLAB. Editing options.json directly has no
% lasting effect -- the next run silently overwrites it. To change what
% gets tested, edit the relevant tests/*.m file.

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

import matlab.unittest.TestSuite
import matlab.unittest.TestRunner
import matlab.unittest.plugins.CodeCoveragePlugin
import matlab.unittest.plugins.codecoverage.CoberturaFormat

% buildIS.m/exploreIS.m/InstanceSpace.m live at the repo root and, unlike
% core/output/utils (self-added to the path the first time InstanceSpace
% is touched -- see InstanceSpace.ensurePathSetup), are never formally on
% the MATLAB search path: they only resolve today via MATLAB's implicit
% current-folder lookup, which is exactly what breaks once
% matlab.unittest shifts the working directory away from the repo root
% during test execution (observed directly in CI -- see tests/testRepoRoot.m
% for the matching fix on the read/write-path side of this same issue).
% This addpath is what actually makes buildIS/InstanceSpace resolvable
% from inside tests/*.m regardless of that shift.
%
% Derived from this file's own location (mfilename), not pwd: this file
% is meant to be run from the repo root, but a caller launching it from
% elsewhere (e.g. `run('/abs/path/test_integration.m')`) would otherwise
% add the wrong directory to the path instead of failing loudly.
repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot);

if ~isfolder('./test/data/')
    mkdir('./test/data/');
end

suite = TestSuite.fromFolder('tests', 'IncludingSubfolders', true);
runner = TestRunner.withTextOutput();

coverageReportFile = 'coverage.xml';
sourceFolders = {'.', 'core', 'output', 'utils'};
runner.addPlugin(CodeCoveragePlugin.forFolder(sourceFolders, ...
    'IncludingSubfolders', false, 'Producing', CoberturaFormat(coverageReportFile)));

results = runner.run(suite);

fprintf('\n[TEST] ================= Summary =================\n');
for i = 1:numel(results)
    if results(i).Passed
        status = 'PASS';
    else
        status = 'FAIL';
    end
    fprintf('[TEST] [%s] %s\n', status, results(i).Name);
end
nPassed = sum([results.Passed]);
nCases  = numel(results);
fprintf('[TEST] %d/%d cases passed.\n', nPassed, nCases);
fprintf('[TEST] Code coverage report written to %s.\n', coverageReportFile);

if nPassed == nCases
    fprintf('EOF:SUCCESS\n');
else
    fprintf('EOF:ERROR\n');
    error('ISA:test_integration:caseFailures', '%d of %d test cases failed.', nCases-nPassed, nCases);
end
