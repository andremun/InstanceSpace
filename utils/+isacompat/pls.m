function [XL,YL,XS,YS,beta,pctvar,mse,stats] = pls(X,Y,components)
% pls  SIMPLS with MATLAB's centered score/weight convention.
% Octave uses the de Jong SIMPLS recurrence; MATLAB retains plsregress.
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

if ~isacompat.isOctave()
    [XL,YL,XS,YS,beta,pctvar,mse,stats] = plsregress(X,Y,components);
    return
end
X0=X-mean(X,1); Y0=Y-mean(Y,1);
if components>rank(X0)
    error('ISA:compat:plsRank','PLS components exceed the rank of centered features.');
end
S=X0'*Y0;
XL=zeros(size(X,2),components); YL=zeros(size(Y,2),components);
XS=zeros(size(X,1),components); W=XL; V=XL; YS=XS;
for k=1:components
    [left,~,~]=svd(S,'econ'); r=left(:,1);
    t=X0*r;
    if k>1
        % Keep scores orthogonal despite accumulated floating-point error.
        correction=XS(:,1:k-1)'*t;
        r=r-W(:,1:k-1)*correction; t=X0*r;
    end
    scale=norm(t);
    if scale<=eps*norm(X0,'fro')
        error('ISA:compat:plsRank','PLS could not form an independent score.');
    end
    r=r/scale; t=t/scale;
    p=X0'*t; q=Y0'*t; v=p;
    if k>1, v=v-V(:,1:k-1)*(V(:,1:k-1)'*v); end
    v=v/norm(v); S=S-v*(v'*S);
    W(:,k)=r; XS(:,k)=t; XL(:,k)=p; YL(:,k)=q; V(:,k)=v;
    YS(:,k)=Y0*q;
end
coeff=W*YL'; beta=[mean(Y,1)-mean(X,1)*coeff;coeff];
pctvar=[sum(XL.^2,1)/sum(X0(:).^2);sum(YL.^2,1)/sum(Y0(:).^2)];
mse=[]; stats=struct('W',W);
end
