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
function X = sobolCandidates(n)
% Two-dimensional Sobol net with a random digital shift in Octave.
% Uses Gray-code direction numbers; intentionally not MATLAB's scramble.
if ~isacompat.isOctave()
    s = scramble(sobolset(2,'Skip',1),'MatousekAffineOwen');
    X = net(s,n); return;
end
v = zeros(2,32,'uint32');
v(1,:) = uint32(2.^(31:-1:0)); v(2,1)=v(1,1);
for j=2:32, v(2,j)=bitxor(v(2,j-1),bitshift(v(2,j-1),-1)); end
shift=uint32(floor(rand(2,1)*2^32)); X=zeros(n,2);
for i=1:n
    g=bitxor(uint32(i),bitshift(uint32(i),-1)); q=shift;
    for j=1:32
        if bitget(g,j), q=bitxor(q,v(:,j)); end
    end
    X(i,:)=double(q)'/2^32;
end
end
