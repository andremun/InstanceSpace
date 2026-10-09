function inside = ISAfootprintContains(geometry, queries, tolerance)
% ISAfootprintContains Closed TRACE membership, shared with pyInstanceSpace.
% Includes the interior and points within the explicit Euclidean boundary
% distance tolerance (projection units). Default zero preserves exact closed
% boundary semantics. The caller owns the error budget; there is no implicit
% epsilon multiplier or dependence on the query batch. Holes remain excluded
% beyond the specified distance. Nonfinite queries are outside.
%
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
% Copyright (c) 2026 Mario Andres Munoz Acosta and contributors
if nargin<3, tolerance=0; end
if ~isnumeric(tolerance) || ~isscalar(tolerance) || ~isreal(tolerance) || ...
        ~isfinite(tolerance) || tolerance<0
    error('ISA:footprint:tolerance','Boundary tolerance must be finite and nonnegative.');
end
inside = false(size(queries,1),1);
if isempty(geometry), return; end
isPoly = isa(geometry,'polyshape');
if isPoly
    dims=2;
    [x,y]=boundary(geometry); vertices=[x y];
    finiteVertices=vertices(all(isfinite(vertices),2),:);
    facets=zeros(0,2);
    starts=[1;find(any(~isfinite(vertices),2))+1];
    stops=[find(any(~isfinite(vertices),2))-1;size(vertices,1)];
    for k=1:numel(starts)
        ids=starts(k):stops(k);
        if numel(ids)>1, facets=[facets;ids(:),ids([2:end 1])']; end %#ok<AGROW>
    end
elseif isacompat.isAlphaShape(geometry)
    [facets,vertices]=boundaryFacets(geometry);
    dims=size(vertices,2);
    finiteVertices=vertices(unique(facets(:)),:);
else
    error('ISA:footprint:geometry','Expected a polyshape or alpha shape.');
end
if ~ismatrix(queries) || size(queries,2)~=dims
    error('ISA:footprint:queries','Query coordinates have the wrong dimensions.');
end
if isempty(finiteVertices) || isempty(facets), return; end
finite=find(all(isfinite(queries),2));
if isempty(finite), return; end
if isPoly
    inside(finite)=isinterior(geometry,queries(finite,:));
elseif isa(geometry,'isacompat.AlphaShape')
    % The compatibility shape's historical barycentric band is not Euclidean.
    inside(finite)=inShape(geometry,queries(finite,:),0);
else
    inside(finite)=inShape(geometry,queries(finite,:));
end
if tolerance==0, return; end
for k=1:size(facets,1)
    facet=vertices(facets(k,:),:);
    selected=find(~inside & all(isfinite(queries),2) & ...
        all(queries>=min(facet,[],1)-tolerance,2) & ...
        all(queries<=max(facet,[],1)+tolerance,2));
    if isempty(selected), continue; end
    q=queries(selected,:); distances=inf(size(q,1),1);
    for j=1:size(facet,1)
        first=facet(j,:); second=facet(mod(j,size(facet,1))+1,:);
        edge=second-first; denominator=sum(edge.^2);
        if denominator==0
            candidate=sqrt(sum((q-first).^2,2));
        else
            fraction=max(0,min(1,(q-first)*edge'/denominator));
            candidate=sqrt(sum((q-first-fraction.*edge).^2,2));
        end
        distances=min(distances,candidate);
    end
    if dims==3
        first=facet(1,:); u=facet(2,:)-first; v=facet(3,:)-first;
        normal=cross(u,v); magnitude=norm(normal);
        if magnitude>0
            normal=normal/magnitude; signed=(q-first)*normal';
            projected=q-first-signed.*normal;
            s=cross(projected,repmat(v,size(q,1),1),2)*normal'/magnitude;
            t=cross(repmat(u,size(q,1),1),projected,2)*normal'/magnitude;
            onTriangle=s>=0 & t>=0 & s+t<=1;
            distances(onTriangle)=min(distances(onTriangle),abs(signed(onTriangle)));
        end
    end
    inside(selected)=distances<=tolerance;
end
end
