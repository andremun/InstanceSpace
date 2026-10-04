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
fprintf('[PORTABLE] PASS: alpha complex area/volume, regions, holes and boundaries.\n');
end
