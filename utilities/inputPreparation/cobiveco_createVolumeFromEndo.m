function [vol,back,shell] = cobiveco_createVolumeFromEndo(endo, numNodesOnHull, lambda, tol, maxit)

if nargin < 2 || isempty(numNodesOnHull)
    numNodesOnHull = 20000;
end
if nargin < 3 || isempty(lambda)
    lambda = 1e0;
end
if nargin < 4 || isempty(tol)
    tol = 1e-8;
end
if nargin < 5 || isempty(maxit)
    maxit = 1000;
end

%%
hull = vtkConvhull(endo);
hull = remeshTriangleMesh(hull, 1000);
hull = vtkConvhull(hull);
hull = remeshTriangleMesh(hull, 1000);
hull = vtkWindowedSincPolyDataFilter(hull, 20);
hull_tr = vtkToTriangulation(hull);
hull.points = hull.points + 15*hull_tr.vertexNormal;
hull = remeshTriangleMesh(hull, numNodesOnHull);
hull = vtkWindowedSincPolyDataFilter(hull, 20);
vtkWrite(hull, 'hull.ply');

back = tetrahedralizeTriangleMesh(hull, 1, 0.7);
meanEdgLen = mean(vtkEdgeLengths(back));

%%
endo = vtkLinearSubdivisionFilter(endo);
ids = vtkMapPointIds(endo, back);
distVec = back.points-endo.points(ids,:);
dist = sqrt(sum((distVec).^2,2));

endo_tr = vtkToTriangulation(endo);
endo_normals = endo_tr.vertexNormal;
dp = dot(distVec, endo_normals(ids,:), 2);
dist = double(sign(dp) .* dist);

if lambda > 0
    L = massmatrix(double(back.points), double(back.cells), 'voronoi') \ cotmatrix(double(back.points), double(back.cells));
    A = speye(size(L))+lambda*(L'*L);
    M = ichol_autocomp(A);
    dist = pcg(A, dist, tol, maxit, M, M', dist);
end

back.pointData.dist = dist;

[layers,mmgStatus1,mmgOutput1] = mmg(back, dist, sprintf('-ls 0 -nr -hausd %1.5e -hmin %1.5e -hmax %1.5e', [0.1 0.9 1.1]*meanEdgLen));
if mmgStatus1 ~= 0
    error(mmgOutput1);
end

%%
layers.pointData.ids = int32(1:size(layers.points,1))';
core = vtkThreshold(layers, 'cells', 'class', [3 3]);
shell = vtkThreshold(layers, 'cells', 'class', [2 2]);
core_sur = vtkDataSetSurfaceFilter(core);
shell_sur = vtkDataSetSurfaceFilter(shell);
ids_core = core_sur.pointData.ids;
ids_shell = shell_sur.pointData.ids;
ids_shell = setdiff(ids_shell, ids_core);
[~,ids1] = ismember(ids_core, shell.pointData.ids);
[~,ids0] = ismember(ids_shell, shell.pointData.ids);

L = cotmatrix(double(shell.points), double(shell.cells));
shell.pointData.laplace = solveLaplace(L, [ids0; ids1], [zeros(size(ids0)); ones(size(ids1))], tol, maxit);

[vol,mmgStatus2,mmgOutput2] = mmg(shell, shell.pointData.laplace, sprintf('-ls 0.5 -nr -hausd %1.5e -hmin %1.5e -hmax %1.5e', [0.1 0.4 0.6]*meanEdgLen));
if mmgStatus2 ~= 0
    error(mmgOutput2);
end

vol = vtkThreshold(vol, 'cells', 'class', [2 2]);

end