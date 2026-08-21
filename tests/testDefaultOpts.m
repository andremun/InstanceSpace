function opts = testDefaultOpts()
% testDefaultOpts  Baseline opts struct shared by every tests/*.m test
% class (#39). Kept close to the published defaults (see README.md /
% docs/tech/isa_refactor_plan_v1.7.pdf Appendix A) except where noted.
% A single shared copy, not one per test class, so the baseline can't
% independently drift between PipelineOptionsTest/ClassApiTest/
% MigrationTest/RegressionTest the way test_integration.m's single
% defaultOpts() never could.

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

opts.general.parallel = false;
opts.general.ncores = 18;

opts.perf.MaxPerf = false;          % True if Y is a performance measure to maximize, False if it is a cost measure to minimise.
opts.perf.AbsPerf = true;           % True if an absolute performance measure, False if a relative performance measure
opts.perf.epsilon = 0.20;           % Threshold of good performance
opts.perf.betaThreshold = 0.55;     % Beta-easy threshold
opts.auto.preproc = true;           % Automatic preprocessing on. Set to false if you don't want any preprocessing
opts.bound.flag = true;             % Bound the outliers. True if you want to bound the outliers, false if you don't
opts.norm.flag = true;              % Normalize/Standarize the data. True if you want to apply Box-Cox and Z transformations to stabilize the variance and scale N(0,1)

opts.selvars.smallscaleflag = false;
opts.selvars.smallscale = 0.50;
opts.selvars.fileidxflag = false;
opts.selvars.fileidx = '';
opts.selvars.densityflag = false;
opts.selvars.mindistance = 0.1;
opts.selvars.type = 'Ftr&Good';

opts.sifted.flag = true;
opts.sifted.rho = 0.1;              % Minimum correlation value acceptable between performance and a feature. Between 0 and 1
% K=5 (not the default 10) so SIFTED's k-means+GA clustering path is
% actually exercised: the reference dataset has only 10 features, and
% K=nfeats triggers SIFTED's "fewer features than clusters" skip branch
% instead (spec §12.1 testing guidance: "SIFTED may be tested with
% opts.sifted.K = 5 or 6").
opts.sifted.K = 5;
opts.sifted.MaxIter = 1000;
opts.sifted.Replicates = 100;

opts.pilot.analytic = false;        % Numerical (BFGS) by default; specific cases below flip this on.
opts.pilot.ntries = 3;              % Small trial count for a quick smoke test; increase for a real run.
opts.pilot.dims = 2;                % 2D by default; specific cases below flip this to 3.
opts.pilot.method = 'standard';     % 'standard' (BFGS/analytic) or 'pls'; specific cases below flip this.
opts.pilot.alpha = 1.0;             % Performance-reconstruction cost weight; specific cases below vary this.

opts.cloister.pval = 0.05;
opts.cloister.corrThreshold = 0.7;

opts.pythia.flag = true;
opts.pythia.classifier = 'knn';      % Registry name: 'knn','svm','tree','nb','linear','ensemble'
opts.pythia.tuning = 'sobol';
opts.pythia.nTuningIter = 5;         % Small evaluation budget for a quick smoke test; increase for a real run.
opts.pythia.kFold = 5;
opts.pythia.ispolykrnl = false;
opts.pythia.useweights = false;
opts.pythia.skip = false;

opts.trace.method = 'trace3';
opts.trace.PI = 0.55;               % Purity threshold

opts.outputs.csv = true;
opts.outputs.web = false;           % NOTE: MAKE THIS FALSE IF YOU ARE USING THIS CODE LOCALLY - This flag is only useful if the system is being used 'online' through matilda.unimelb.edu.au
opts.outputs.png = true;
end
