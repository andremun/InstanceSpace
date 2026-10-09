function runBoundarySmoke()
% Explicit distance semantics for the toolkit alpha shape on MATLAB/Octave.
for dims=[2 3]
    vertices=[zeros(1,dims);eye(dims)];
    shape=isacompat.AlphaShape(vertices); shape.Alpha=Inf;
    q=zeros(4,dims); q(:,1)=0.25;
    q(1,end)=0; q(2,end)=-1/2048; q(3,end)=-1/512;
    q(4,:)=-3/4096;
    assert(isequal(ISAfootprintContains(shape,q),logical([1;0;0;0])));
    assert(isequal(ISAfootprintContains(shape,q,1/1024),logical([1;1;0;0])));
end
fprintf('[PORTABLE] PASS: explicit boundary distance in 2D and 3D.\n');
end
