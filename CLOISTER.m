function out = CLOISTER(X, A, opts)
% -------------------------------------------------------------------------
% CLOISTER.m
% -------------------------------------------------------------------------
%
% By: Mario Andres Munoz Acosta
%     School of Mathematics and Statistics
%     The University of Melbourne
%     Australia
%     2020
%
% -------------------------------------------------------------------------

fprintf('  -> CLOISTER is using correlation to estimate a boundary for the space.\n');

nfeats = size(X,2);
[rho,pval] = corr(X);
rho = rho.*(pval<opts.pval);

Xbnds = [min(X); max(X)];
% Guard: if too many features, de2bi would produce an intractable matrix.
% Use convex hull of Z as a safe fallback.
MAX_FEATS = 20;
if nfeats > MAX_FEATS
    warning('ISA:CLOISTER:tooManyFeatures', ...
        'CLOISTER skipped: %d features exceeds limit of %d. Using convex hull as boundary.', ...
        nfeats, MAX_FEATS);
    Zall = X*A';
    Kedge = convhull(Zall(:,1), Zall(:,2));
    out.Zedge  = Zall(Kedge,:);
    out.Zecorr = out.Zedge;
    fprintf('  -> CLOISTER has completed.\n');
    return;
end
% Pure-MATLAB replacement for de2bi (no Communications Toolbox required)
idx = rem(floor((0:2^nfeats-1)' .* 2.^(-(nfeats-1:-1:0))), 2) + 1;
ncomb = size(idx,1);
Xedge = zeros(ncomb,nfeats);
remove = false(ncomb,1);
for i=1:ncomb
   ind = sub2ind([2 nfeats],idx(i,:),1:nfeats);
   Xedge(i,:) = Xbnds(ind)';
   for j=1:nfeats
       for k=j+1:nfeats
           % Check for valid points give the correlation trend
           if rho(j,k)>opts.cthres && sign(Xedge(i,j))~=sign(Xedge(i,k))
               remove(i) = true;
           elseif rho(j,k)<-opts.cthres && sign(Xedge(i,j))==sign(Xedge(i,k))
               remove(i) = true;
           end
           if remove(i)
               break;
           end
       end
       if remove(i)
           break;
       end
   end
end
Zedge = Xedge*A';
Kedge = convhull(Zedge(:,1),Zedge(:,2));
out.Zedge = Zedge(Kedge,:);

try
    Xecorr = Xedge(~remove,:);
    Zecorr = Xecorr*A';
    Kecorr = convhull(Zecorr(:,1),Zecorr(:,2));
    out.Zecorr = Zecorr(Kecorr,:);
catch
    fprintf('  -> The acceptable correlation threshold was too strict.\n');
    fprintf('  -> The features are weakely correlated.\n');
    fprintf('  -> Please consider increasing it.\n');
    out.Zecorr = out.Zedge;
end
fprintf('-------------------------------------------------------------------------\n');
fprintf('  -> CLOISTER has completed.\n');