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
function runOptimizationSmoke()
t=(1:30)'; X=[sin(t),cos(t),t/30,sin(t/3)]; Y=X*[1 2;2 -1;3 1;-1 2];
labels={'a','b','c','d'};
for dims=[2 3]
    opts=struct('dims',dims,'method','pls','verbose',false,'ntries',1,'seed',9);
    state=rng; out=PILOT(X,Y,labels,opts); assert(isequal(rng,state));
    assert(norm(out.Z-(X-mean(X))*out.A','fro')<1e-10);
    assert(norm(out.Z'*out.Z-eye(dims),'fro')<1e-10);
    opts.method='standard'; opts.analytic=false;
    out=PILOT(X,Y,labels,opts); assert(isequal(rng,state));
    assert(all(isfinite(out.Z(:))) && isfinite(out.error));
    assert(norm(out.Z-X*out.A','fro')<1e-10);
end
if isacompat.isOctave()
    settings=struct('PopulationSize',50,'MaxGenerations',10,'FitnessLimit',0, ...
        'FunctionTolerance',1e-3,'MaxStallGenerations',3);
    [best,info]=isacompat.integerSearch(@(x)sum((x-[2 3]).^2),[3 4],settings);
    assert(isequal(best,[2 3]) && info.fitness==0);
    rng(9); [best,info]=isacompat.integerSearch(@(x)sum((x-1).^2),[10 10],settings);
    assert(all(best>=1 & best<=10 & best==round(best)) && isfinite(info.fitness));
    rng(9); A=isacompat.sobolCandidates(15); rng(9); B=isacompat.sobolCandidates(15);
    assert(isequal(A,B) && all(A(:)>=0 & A(:)<1));
    assert(size(unique(floor(A*4),'rows'),1)==15);
end
fprintf('[PORTABLE] PASS: numerical/PLS projections and bounded integer/Sobol search.\n');
end
