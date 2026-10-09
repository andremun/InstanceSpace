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
classdef AlphaShape
% Regularized full-dimensional Delaunay alpha complex for TRACE.
% Regions connect through facets. Holes are retained; boundary is inclusive.
properties (SetAccess=private)
    Points
    Simplices
    Radii
    Measures
end
properties (Access=private)
    SpatialIndex=[]
end
properties
    Alpha=Inf
    RegionThreshold=0
end
methods
    function obj=AlphaShape(P)
        obj.Points=unique(P,'rows'); P=obj.Points; d=size(P,2);
        if ~ismember(d,[2 3]) || any(~isfinite(P(:)))
            error('ISA:compat:geometry','Expected finite 2D or 3D points.');
        end
        obj.Simplices=zeros(0,d+1); obj.Radii=[]; obj.Measures=[];
        if size(P,1)<=d || rank(P-P(1,:))<d, return; end
        T=delaunayn(P,{'Qt','Qbb','Qc','Qz'}); r=zeros(size(T,1),1); m=r;
        for k=1:size(T,1)
            E=P(T(k,2:end),:)-P(T(k,1),:);
            m(k)=abs(det(E))/factorial(d);
            if rcond(E)>eps, r(k)=norm(E\(sum(E.^2,2)/2)); else, r(k)=Inf; end
        end
        keep=m>0 & isfinite(r); obj.Simplices=T(keep,:);
        obj.Radii=r(keep); obj.Measures=m(keep);
        lower=Inf(size(obj.Simplices,1),d); upper=-lower;
        for j=1:d+1
            vertices=P(obj.Simplices(:,j),:);
            lower=min(lower,vertices); upper=max(upper,vertices);
        end
        obj.SpatialIndex=buildIndex(lower,upper,(1:size(lower,1))');
        nearest=Inf(size(P,1),1);
        for k=1:numel(obj.Radii)
            ids=obj.Simplices(k,:); nearest(ids)=min(nearest(ids),obj.Radii(k));
        end
        obj.Alpha=max(nearest);
    end
    function a=alphaSpectrum(obj), a=sort(unique(obj.Radii),'descend'); end
    function [T,regions,measure,ids]=active(obj)
        ids=find(obj.Radii<=obj.Alpha*(1+64*eps)); T=obj.Simplices(ids,:);
        measure=obj.Measures(ids); n=size(T,1); regions=zeros(n,1);
        if n==0, return; end
        d=size(T,2)-1; F=[]; owners=[];
        for j=1:d+1
            F=[F;sort(T(:,setdiff(1:d+1,j)),2)]; owners=[owners;(1:n)'];
        end
        [~,~,groups]=unique(F,'rows'); [sorted,order]=sort(groups);
        pair=find(diff(sorted)==0); a=owners(order(pair)); b=owners(order(pair+1));
        adj=sparse([a;b],[b;a],1,n,n); count=0;
        for k=1:n
            if regions(k), continue; end
            count=count+1; queue=k; regions(k)=count; pos=1;
            while pos<=numel(queue)
                next=find(adj(queue(pos),:)); next=next(regions(next)==0);
                regions(next)=count; queue=[queue,next]; pos=pos+1;
            end
        end
        totals=accumarray(regions,measure); [~,ord]=sort(totals,'descend');
        mapping=zeros(count,1); valid=ord(totals(ord)>=obj.RegionThreshold);
        mapping(valid)=1:numel(valid); regions=mapping(regions);
        keep=regions>0; T=T(keep,:); regions=regions(keep); measure=measure(keep); ids=ids(keep);
    end
    function n=numRegions(obj), [~,r]=obj.active(); n=max([0;r]); end
    function a=area(obj), [~,~,m]=obj.active(); a=sum(m); end
    function a=volume(obj), a=area(obj); end
    function [F,P]=boundaryFacets(obj,region)
        [T,r]=obj.active(); P=obj.Points;
        if nargin>1, T=T(r==region,:); end
        d=size(P,2); F=zeros(0,d);
        for j=1:d+1, F=[F;sort(T(:,setdiff(1:d+1,j)),2)]; end
        if isempty(F), return; end
        [U,~,g]=unique(F,'rows'); counts=accumarray(g,1); F=U(counts==1,:);
    end
    function inside=inShape(obj,Q,tolerance)
        if nargin<3, tolerance=1e-10; end
        [~,~,~,activeIds]=obj.active(); inside=false(size(Q,1),1);
        if isempty(activeIds) || isempty(Q), return; end
        enabled=false(size(obj.Simplices,1),1); enabled(activeIds)=true;
        inside=queryIndex(obj.SpatialIndex,obj.Points,obj.Simplices,enabled,Q, ...
            (1:size(Q,1))',inside,tolerance);
    end
    function h=plot(obj,varargin)
        if size(obj.Points,2)==2, F=obj.active(); else, F=boundaryFacets(obj); end
        h=patch('Faces',F,'Vertices',obj.Points,varargin{:});
    end
end
end

function node=buildIndex(lower,upper,ids)
% A reusable bounding-volume hierarchy over the immutable Delaunay complex.
node=struct('lower',min(lower(ids,:),[],1),'upper',max(upper(ids,:),[],1), ...
    'ids',[],'left',[],'right',[]);
if numel(ids)<=16
    node.ids=ids; return;
end
centers=lower(ids,:)/2+upper(ids,:)/2;
[~,axis]=max(max(centers,[],1)-min(centers,[],1));
[~,order]=sort(centers(:,axis)); middle=floor(numel(ids)/2);
node.left=buildIndex(lower,upper,ids(order(1:middle)));
node.right=buildIndex(lower,upper,ids(order(middle+1:end)));
end
function inside=queryIndex(node,P,T,enabled,Q,ids,inside,tol)
ids=ids(~inside(ids) & all(Q(ids,:)>=node.lower-tol & Q(ids,:)<=node.upper+tol,2));
if isempty(ids), return; end
if isempty(node.ids)
    inside=queryIndex(node.left,P,T,enabled,Q,ids,inside,tol);
    inside=queryIndex(node.right,P,T,enabled,Q,ids,inside,tol);
else
    for k=node.ids(enabled(node.ids))'
        V=P(T(k,:),:);
        candidates=ids(~inside(ids) & all(Q(ids,:)>=min(V)-tol & Q(ids,:)<=max(V)+tol,2));
        if isempty(candidates), continue; end
        B=(Q(candidates,:)-V(1,:))/(V(2:end,:)-V(1,:));
        inside(candidates)=all(B>=-tol,2) & sum(B,2)<=1+tol;
    end
end
end
