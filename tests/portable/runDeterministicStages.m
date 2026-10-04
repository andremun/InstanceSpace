function runDeterministicStages()
% runDeterministicStages  Portable preprocessing and hull contracts.
% Called by runPortableSmoke after it sets up the core/utils search paths.
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

X = [1 5 NaN 3; 2 NaN NaN 3; 3 NaN NaN NaN; 4 NaN NaN 3];
[med, spread] = isacompat.columnQuartiles(X);
assert(isequaln(med, [2.5 5 NaN 3]));
assert(isequaln(spread, [2 0 NaN 0]));
[lo, hi] = isacompat.columnExtrema(X);
assert(isequaln(lo, [1 5 NaN 3]));
assert(isequaln(hi, [4 5 NaN 3]));
[med, spread] = isacompat.columnQuartiles([2 NaN -3]);
assert(isequaln(med, [2 NaN -3]));
assert(isequaln(spread, [0 NaN 0]));
[lo, hi] = isacompat.columnExtrema([-Inf NaN; Inf NaN]);
assert(isequaln(lo, [-Inf NaN]) && isequaln(hi, [Inf NaN]));

% For n=4 the exact two-sided Pearson significance is 1-abs(r).
u = [-1;-1;1;1]/2;
v = [-1;1;-1;1]/2;
y = 0.5*u + sqrt(0.75)*v;
[r,p] = isacompat.pearsonColumns([u y -y v ones(4,1) [NaN;1;2;3]]);
assert(abs(r(1,2)-0.5)<1e-14 && abs(r(1,3)+0.5)<1e-14);
assert(abs(p(1,2)-0.5)<1e-14 && abs(p(1,3)-0.5)<1e-14);
assert(abs(p(1,4)-1)<1e-14);
assert(all(all(isnan(r(5:6,:)))) && all(all(isnan(p(:,5:6)))));
[r,p] = isacompat.pearsonColumns([u -u]);
assert(abs(r(1,2)+1)<1e-14 && p(1,2)<1e-7);

opts = struct('MaxPerf',false,'AbsPerf',true,'epsilon',0.5, ...
    'betaThreshold',0.5,'auto',true,'bound',true,'norm',false,'iqrMultiplier',1);
X = [0 1; 1 2; 2 NaN; 3 4; 100 5];
Y = [0.1 0.9; 0.6 0.2; 0.4 0.8; 0.8 0.7; NaN NaN];
[Xn,Yn,trained] = PRELIM(X,Y,opts);
assert(isequaln(Yn,Y));
assert(isequal(trained.Ybin, logical([1 0;0 1;1 0;0 0;0 0])));
assert(isequal(trained.P,[1;2;1;2;0]));
assert(isequaln(trained.Ybest,[0.1;0.2;0.4;0.7;NaN]));
assert(isequal(trained.numGoodAlgos,[1;1;1;0;0]));
assert(~any(trained.beta));
assert(isequaln(trained.medval,[2 3]));
assert(isequaln(trained.iqrange,[26.5 3]));
assert(Xn(5,1)==28.5 && isnan(Xn(3,2)));
% Applying the fitted bounds reproduces the training transform.
oldWarning = warning('off','ISA:InstanceSpace:outOfDistribution');
cleanup = onCleanup(@() warning(oldWarning)); %#ok<NASGU>
[Xe,Ye] = PRELIM(X,Y,opts,trained);
assert(isequaln(Xe,Xn) && isequaln(Ye,Yn));
% Exercise minimization/maximization with absolute and relative thresholds.
rawPerf = [1 2;4 2;3 5;NaN 1];
modeOpts = opts; modeOpts.auto = false;
expectedMasks = {logical([1 1;0 1;0 0;0 1]), ...
    logical([0 1;1 1;1 1;0 0]), logical([1 0;0 1;1 0;0 1]), ...
    logical([1 1;1 1;1 1;0 1])};
for mode = 1:4
    modeOpts.MaxPerf = mod(mode,2)==0;
    modeOpts.AbsPerf = mode<=2;
    if modeOpts.AbsPerf, modeOpts.epsilon = 2; else, modeOpts.epsilon = 0.5; end
    [~,transformed,labels] = PRELIM((1:4)',rawPerf,modeOpts);
    assert(isequal(labels.Ybin,expectedMasks{mode}));
    assert(isequal(labels.beta,sum(expectedMasks{mode},2)>1));
    if mode<=2
        assert(isequaln(transformed,rawPerf));
    elseif mode==3
        assert(isequaln(transformed,[0 1;1 0;0 5/3-1;NaN 0]));
    else
        assert(isequaln(transformed,[0.5 0;0 0.5;1-3/5 0;NaN 0]));
    end
end
% The evaluation-only Box-Cox path accepts known fitted parameters.
fit = struct('lambdaX',[0 1],'minX',[0 0],'muX',[0 0], ...
    'sigmaX',[1 2],'lambdaY',[0 1],'minY',0,'muY',[0 0],'sigmaY',[1 2]);
evalOpts = opts; evalOpts.norm = true; evalOpts.bound = false;
rawX = [0 0;1 1;3 3]; rawY = [1 2;2 3;4 5];
[Xe,Ye] = PRELIM(rawX,rawY,evalOpts,fit);
assert(norm(Xe-[log(rawX(:,1)+1),rawX(:,2)/2],'fro')<1e-14);
expectedY = [log(rawY(:,1)+eps),(rawY(:,2)+eps-1)/2];
assert(norm(Ye-expectedY,'fro')<1e-14);
if isacompat.isOctave()
    [fittedX,~,fitted] = PRELIM(rawX,rawY,evalOpts);
    appliedX = PRELIM(rawX,rawY,evalOpts,fitted);
    assert(norm(fittedX-appliedX,'fro')<1e-10);
    assertStageError(@() isacompat.pearsonColumns([1 2;3 4]), ...
        'ISA:compat:insufficientCorrelationRows');
end

% Exact rectangle/cuboid geometry, without depending on vertex order.
[a,b,c] = ndgrid([-1 1],[-2 2],[-3 3]);
X = [a(:), b(:), c(:)];
opts = struct('pval',0.05,'corrThreshold',0.7,'maxFeatures',20);
out = CLOISTER(X,[1 0 0;0 1 0],opts);
assert(abs(polyarea(out.Zedge(:,1),out.Zedge(:,2))-8)<1e-12);
assert(isequal(out.Zedge(1,:),out.Zedge(end,:)));
assert(isequal(out.Zedge,out.Zecorr));
assert(isempty(out.ZedgeFaces));
out = CLOISTER(X,eye(3),opts);
[~,volume] = convhull(out.Zedge);
assert(abs(volume-48)<1e-10);
assert(size(out.ZedgeFaces,2)==3);
assert(all(out.ZedgeFaces(:)>=1 & out.ZedgeFaces(:)<=size(out.Zedge,1)));
assert(norm(mean(out.Zedge,1))<1e-12);
% PLS centering translates the hull; it must not change the volume.
shift = [2 4 8];
shifted = CLOISTER(X,eye(3),opts,shift);
assert(norm(mean(shifted.Zedge,1)+shift)<1e-12);
% Coplanar 3D projection returns a triangulated flat polygon.
A = [1 0;0 1;1 1];
flat = CLOISTER(X(:,1:2),A,opts);
assert(max(abs(flat.Zedge(:,3)-flat.Zedge(:,1)-flat.Zedge(:,2)))<1e-12);
assert(size(flat.ZedgeFaces,1)==size(flat.Zedge,1)-2);
% Correlation restrictions collapse the valid corners to a line: use full hull.
x = [-2;-1;1;2];
restricted = CLOISTER([x 2*x],eye(2),opts);
assert(isequal(restricted.Zedge,restricted.Zecorr));
assert(abs(polyarea(restricted.Zedge(:,1),restricted.Zedge(:,2))-32)<1e-12);
% A degenerate full projection retains the stage's diagnostic error.
assertStageError(@() CLOISTER(X,[1 1 1;2 2 2],opts), ...
    'ISA:CLOISTER:degenerateBoundary');
% Missing values must not contaminate feature extrema.
missingX = X; missingX(1,1) = NaN;
missing = CLOISTER(missingX,[1 0 0;0 1 0],opts);
assert(abs(polyarea(missing.Zedge(:,1),missing.Zedge(:,2))-8)<1e-12);
% The feature-count fallback encloses projected observations.
opts.maxFeatures = 2;
oldLimit = warning('off','ISA:CLOISTER:tooManyFeatures');
limitCleanup = onCleanup(@() warning(oldLimit)); %#ok<NASGU>
fallback = CLOISTER(X,eye(3),opts);
[~,volume] = convhull(fallback.Zedge);
assert(abs(volume-48)<1e-10 && isequal(fallback.Zedge,fallback.Zecorr));
fprintf('[PORTABLE] PASS: PRELIM bounds/labels/evaluation and CLOISTER 2D/3D.\n');
end

function assertStageError(fcn, identifier)
try
    fcn();
catch err
    assert(strcmp(err.identifier,identifier), ...
        'Expected %s; received %s: %s',identifier,err.identifier,err.message);
    return
end
error('ISA:portable:missingError','Expected error %s.',identifier);
end
