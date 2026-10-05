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
function runGeometrySmoke()
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
oldPath=path; cleanup=onCleanup(@() path(oldPath)); %#ok<NASGU>
addpath(fullfile(root,'utils'));
state=rng; rngCleanup=onCleanup(@() rng(state)); %#ok<NASGU>
rng(73);
for dims=[2 3]
    P=dec2bin(0:2^dims-1)-'0'; shape=isacompat.AlphaShape(P);
    assert(abs(area(shape)-1)<1e-12);
    assert(all(inShape(shape,[P;ones(1,dims)/2])));
    assert(~inShape(shape,2*ones(1,dims)));
    assert(numRegions(shape)==1);
    small=isacompat.AlphaShape([zeros(1,dims);eye(dims)]);
    assert(abs(area(small)-1/factorial(dims))<1e-12);
    joined=isacompat.AlphaShape([P;P+4]); joined.Alpha=sqrt(dims);
    assert(numRegions(joined)==2 && abs(area(joined)-2)<1e-12);
    joined.RegionThreshold=1.1; assert(area(joined)==0);
end
shape=isacompat.AlphaShape([0 0;1 0;2 0]); assert(area(shape)==0);
shape=isacompat.AlphaShape([0 0;1 0;1 1;0 1;0 0]); assert(area(shape)==1);
% An annular point set must not fill its central hole at small alpha.
t=(0:15)'*2*pi/16; P=[cos(t),sin(t);2*cos(t),2*sin(t)];
shape=isacompat.AlphaShape(P); shape.Alpha=.8;
assert(area(shape)>0 && ~inShape(shape,[0 0]));
checkIndexedContainment();
fprintf('[PORTABLE] PASS: alpha complex area/volume, regions, holes and boundaries.\n');
end

function checkIndexedContainment()
for dims=[2 3]
    shape=isacompat.AlphaShape(rand(80,dims));
    P=shape.Points; T=shape.Simplices;
    % Vertices, facet centroids and small perturbations test inclusive boundaries.
    faces=zeros(size(T,1),dims);
    for j=1:dims, faces=faces+P(T(:,j),:)/dims; end
    Q=[rand(300,dims)*1.4-.2;P;faces;faces+1e-11;faces-1e-11;NaN(1,dims)];
    radii=sort(shape.Radii);
    for fraction=[.2 .6 1]
        shape.Alpha=radii(max(1,ceil(fraction*numel(radii))));
        for threshold=[0 .001 Inf]
            shape.RegionThreshold=threshold;
            assert(isequal(inShape(shape,Q),exhaustiveContainment(shape,Q)));
        end
    end
    assert(isempty(inShape(shape,zeros(0,dims))));
end
% Representative timing, reported rather than asserted on shared CI hardware.
shape=isacompat.AlphaShape(rand(500,3)); Q=rand(2000,3);
tic; expected=exhaustiveContainment(shape,Q); referenceTime=toc;
tic; actual=inShape(shape,Q); indexedTime=toc;
assert(isequal(actual,expected));
fprintf('[PORTABLE] Containment benchmark: exhaustive %.3fs, indexed %.3fs (%d tetrahedra).\n', ...
    referenceTime,indexedTime,size(shape.Simplices,1));
end
function inside=exhaustiveContainment(shape,Q)
% Independent reference retaining the original exhaustive containment contract.
T=shape.active(); P=shape.Points; inside=false(size(Q,1),1); tol=1e-10;
for k=1:size(T,1)
    V=P(T(k,:),:); ids=find(~inside & all(Q>=min(V)-tol & Q<=max(V)+tol,2));
    if isempty(ids), continue; end
    B=(Q(ids,:)-V(1,:))/(V(2:end,:)-V(1,:));
    inside(ids)=all(B>=-tol,2) & sum(B,2)<=1+tol;
end
end
