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

disp('  -> CLOISTER is using correlation to estimate a boundary for the space.');

nfeats = size(X,2);
[rho,pval] = corr(X);
rho = rho.*(pval<opts.pval);

Xbnds = [min(X); max(X)];
% if no feature selection. then make a note in the boundary construction
% that it won't work, because nfeats is soo large that de2bi wont be able
% to make a matrix.
idx = de2bi(0:2^nfeats-1)+1;
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
    disp('  -> The acceptable correlation threshold was too strict.');
    disp('  -> The features are weakely correlated.')
    disp('  -> Please consider increasing it.');
    out.Zecorr = out.Zedge;
end
disp('-------------------------------------------------------------------------');
disp('  -> CLOISTER has completed.');

