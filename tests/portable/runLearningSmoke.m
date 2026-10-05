function runLearningSmoke()
% runLearningSmoke  Shared untuned KNN contracts; requires Statistics APIs.
% In Octave first run setupOctaveValidation. MATLAB uses installed toolboxes.
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
addpath(fullfile(root,'utils'),fullfile(root,'core'));
if isempty(which('fitcknn')) || isempty(which('cvpartition'))
    error('ISA:compat:missingPackages','KNN tests require Statistics APIs; in Octave run setupOctaveValidation first.');
end
state = rng;
rngCleanup = onCleanup(@() rng(state)); %#ok<NASGU>
t = (1:12)'/12;
% Constant second coordinate also exercises the fitted unit-scale fallback.
Z = [-5-t zeros(size(t));5+t zeros(size(t))];
good = Z(:,1)>0;
Y = [1+3*double(~good),1+2*double(good),ones(size(good))];
Ybin = [good ~good true(size(good))];
labels = {'right','left','constant'};
opts = struct('classifier','knn','tuning','none','params',[3 1;3 1;1 1], ...
    'kFold',3,'verbose',false,'parallel',false,'seed',17);
for weighted = [false true]
    opts.useweights = weighted;
    before = rng;
    out = PYTHIA(Z,Y,Ybin,min(Y,[],2),labels,opts);
    assert(isequal(rng,before));
    assert(isequal(out.Yhat,Ybin) && isequal(out.Ysub,Ybin));
    assert(all(out.accuracy==1) && all(out.precision==1));
    assert(max(abs(out.Pr0hat(:)-double(~Ybin(:))))<1e-12);
    assert(max(abs(out.Pr0sub(:)-double(~Ybin(:))))<1e-12);
    assert(all(out.Pr0subIsProbability(:)));
    assert(all(strcmp(out.scoreType,'probability')));
    assert(out.classifiers{3}.constant && isempty(out.cp{3}));
    for algorithm = 1:2
        cp = out.cp{algorithm};
        coverage = zeros(size(good));
        for fold = 1:cp.NumTestSets
            trainMask = training(cp,fold); testMask = test(cp,fold);
            assert(isequal(size(trainMask),size(good)));
            assert(all(xor(trainMask,testMask)));
            assert(any(good(trainMask)) && any(~good(trainMask)));
            coverage = coverage+testMask;
        end
        assert(all(coverage==1));
        assert(isequal(logical(out.classifiers{algorithm}.ClassNames(:)),[false;true]));
    end
    if weighted
        expectedWeights = abs(Y-min(Y,[],2));
        expectedWeights(expectedWeights==0)=2;
        assert(isequal(out.W,expectedWeights));
    end
    % Held-out points, then missing outcomes: predictions must be preserved.
    Ztest = [-5.5 0;5.5 0];
    expected = logical([0 1 1;1 0 1]);
    outcomes = [4 1 1;1 3 1];
    evaluated = PYTHIA(Ztest,outcomes,expected,[1;1],labels,opts,out);
    assert(isequal(evaluated.Yhat,expected));
    assert(max(abs(evaluated.Pr0hat(:)-double(~expected(:))))<1e-12);
    outcomes(:,1)=NaN; missingLabels=expected; missingLabels(:,1)=false;
    missing = PYTHIA(Ztest,outcomes,missingLabels,[1;1],labels,opts,out);
    assert(isequal(missing.Yhat,expected));
    assert(isnan(missing.accuracy(1)) && missing.accuracy(2)==1);
end
opts.tuning='sobol'; opts.nTuningIter=4; opts.params=[];
before=rng; tuned=PYTHIA(Z,Y,Ybin,min(Y,[],2),labels,opts);
assert(isequal(rng,before));
assert(all(isfinite(tuned.accuracy)));
assert(isequal(size(tuned.tuningCandidates{1}),[4 2]));
assert(all(tuned.Pr0hat(:)>=0 & tuned.Pr0hat(:)<=1));
if isacompat.isOctave()
    unsupported = opts; unsupported.tuning = 'bayes';
    assertLearningError(@() PYTHIA(Z,Y,Ybin,min(Y,[],2),labels,unsupported), ...
        'ISA:compat:unsupportedFeature');
    unsupported = opts; unsupported.classifier = 'svm';
    assertLearningError(@() PYTHIA(Z,Y,Ybin,min(Y,[],2),labels,unsupported), ...
        'ISA:compat:unsupportedFeature');
end
fprintf('[PORTABLE] PASS: weighted/unweighted PYTHIA KNN, CV, none/Sobol tuning and evaluation.\n');
end

function assertLearningError(fcn, identifier)
try
    fcn();
catch err
    assert(strcmp(err.identifier,identifier), ...
        'Expected %s; received %s: %s',identifier,err.identifier,err.message);
    return
end
error('ISA:portable:missingError','Expected error %s.',identifier);
end
