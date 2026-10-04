function report = runPortableSmoke()
% runPortableSmoke  Shared MATLAB/Octave checks without matlab.unittest.
% Run from the repository root: addpath('tests/portable'); runPortableSmoke
% Uses an independent SVD reconstruction oracle, not historical MATLAB data.
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

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
oldPath = path;
cleanup = onCleanup(@() path(oldPath)); %#ok<NASGU>
addpath(fullfile(root, 'utils'), fullfile(root, 'core'));
report = isacompat.diagnostics();
assert(strcmp(report.version, version));
assert(report.apis.corr.available);

% No random input: reproducibility does not depend on engine RNG streams.
t = (1:24)';
% Keep every column nonconstant: correlation of a numerically almost
% constant QR column is ill-conditioned and unsuitable for a strict R2 test.
[Q, ~] = qr([t, t.^2, sin(t), cos(t)], 0);
X = Q * diag([1 2 3 4]);
Y = X * [2 -1; 1 3; -2 2; 3 1];
labels = {'one','two','three','four'};
uniformState = rand('state');
normalState = randn('state');
for dims = [2 3]
    for weight = [0.5 1 3]
        opts = struct('analytic', true, 'method', 'standard', ...
            'parallel', false, 'verbose', false, 'dims', dims, 'alpha', weight);
        out = PILOT(X, Y, labels, opts);
        weighted = [X sqrt(weight)*Y];
        [U, S, V] = svd(weighted, 'econ');
        expected = U(:,1:dims)*S(1:dims,1:dims)*V(:,1:dims)';
        expected(:,5:end) = expected(:,5:end) / sqrt(weight);
        actual = out.Z * [out.B; out.C']';
        assert(isequal(size(out.A), [dims 4]));
        assert(isequal(size(out.Z), [24 dims]));
        assert(norm(out.Z-X*out.A', 'fro') < 1e-11);
        assert(norm(actual-expected, 'fro') < 1e-10*max(1,norm(expected,'fro')));
        expectedError = sum(sum(([X Y]-expected).^2));
        assert(abs(out.error-expectedError) < 1e-10*max(1,expectedError));
        for k = 1:6
            combined = [X Y];
            r = corrcoef(combined(:,k), expected(:,k));
            assert(abs(out.R2(k)-r(1,2)^2) < 1e-10);
        end
        assert(isequal(out.summary(1,2:end), labels));
        displayed = cell2mat(out.summary(2:end,2:end));
        assert(max(abs(displayed(:)-out.A(:))) <= 0.5e-4+eps);
    end
end
assert(isequal(rand('state'), uniformState));
assert(isequal(randn('state'), normalState));
assert(isequal(isacompat.roundDecimal([-1.25 0 1.25],1),[-1.3 0 1.3]));
large = isacompat.roundDecimal([realmax Inf -Inf NaN],4);
assert(large(1)==realmax && large(2)==Inf && large(3)==-Inf && isnan(large(4)));
opts = struct('analytic', true, 'verbose', false);
assertError(@() PILOT(X,Y,labels,struct('dims',4)), 'ISA:PILOT:invalidDims');
assertError(@() PILOT(X,Y,labels,struct('method','invalid')), 'ISA:PILOT:invalidMethod');
badX = X; badX(1) = NaN;
assertError(@() PILOT(badX,Y,labels,opts), 'ISA:PILOT:incompleteData');
if report.isOctave
    assertError(@() PILOT(X,Y,labels,struct('analytic',true,'parallel',true)), ...
        'ISA:compat:unsupportedFeature');

end
runDeterministicStages();
fprintf('[PORTABLE] PASS: %s %s; analytic PILOT 2D/3D, weights, summaries and errors.\n', ...
    report.engine, report.version);
end

function assertError(fcn, identifier)
try
    fcn();
catch err
    assert(strcmp(err.identifier, identifier), ...
        'Expected %s; received %s: %s', identifier, err.identifier, err.message);
    return
end
error('ISA:portable:missingError', 'Expected error %s.', identifier);
end
