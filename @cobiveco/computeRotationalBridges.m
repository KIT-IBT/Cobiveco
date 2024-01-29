function computeRotationalBridges(o)

% Computes the anterior-posterior coordinate on the intervalvular regions
% of the given object,
% i.e. the trajectory distance from posterior to anterior boundary surface.
% The coordinate values are added in the fields "rt", "rtSin" and "rtCos" to the output objects.
%
% computeRotationalBridges(o)
%
% Input:
%   o: instance of the class cobiveco
%
% Output:
%   o.m0RvBridge, object of class cobiveco [struct] (For details see cobiveco class documentation)
%   o.m0LvBridge, object of class cobiveco [struct] (For details see cobiveco class documentation)


if ~o.available.creategraphandsplitmesh
    o.createGraphandSplitMesh;
end

o.printStatus('Computing anterior posterior coordinate in intervalvular regions.');
t = toc;

indicesMesh = knnsearch(o.m0Ventricles.sur.points,o.m0RvBridge.sur.points,  'NSMethod','kdtree');
o.m0Ventricles.sur.pointData.remove = repmat(uint8(0),size(o.m0Ventricles.sur.points,1),1);
indicesMesh = unique(indicesMesh);
o.m0Ventricles.sur.pointData.remove(indicesMesh) = 1;
ventricularmeshBoundary  = vtkDeleteDataArrays(vtkThreshold(o.m0Ventricles.sur, 'points', 'remove', [1 1]));
vtkWrite(ventricularmeshBoundary, sprintf('%sventricularmeshSurf.vtp', o.cfg.outPrefix));

% Create FreeSeptalWall boundary surfaces
[o.m0RvBridge,o.rvBridgeBoundaryFw,o.rvBridgeBoundarySeptum] = createBoundarySurfacesAB(o.m0RvBridge, o.m0Ventricles, o.leftRightAx, -o.leftRightAx, o.longAx, o.cfg.outPrefix);
[o.m0LvBridge,o.lvBridgeBoundaryFw,o.lvBridgeBoundarySeptum] = createBoundarySurfacesAB(o.m0LvBridge, o.m0Ventricles, -o.leftRightAx, o.leftRightAx, o.longAx, o.cfg.outPrefix);

% Create AnteriorPosterior boundary surfaces
[rvOTBoundary, rvITBoundary] = createBoundarySurfacesRT(o.m0RvBridge, o.surPv, o.surTv, o.cfg.outPrefix);
[lvOTBoundary, lvITBoundary] = createBoundarySurfacesRT(o.m0LvBridge, o.surAv, o.surMv, o.cfg.outPrefix);

% create new base boundaries
[o.m0Ventricles, o.ventricularmeshRvBoundary,o.ventricularmeshLvBoundary] = createBoundarySurfaceRVLVOnVentricle(o.m0Ventricles, o.rvBridgeBoundaryFw, o.rvBridgeBoundarySeptum,o.lvBridgeBoundaryFw,o.lvBridgeBoundarySeptum, o.cfg.outPrefix);

% map tv and tm onto bridges and ventricles
o.m0Ventricles = mapCoordinateFromMeshToSubmeshBary(o.m0,o.m0Ventricles, "tv", "tm");
o.m0LvBridge = mapCoordinateFromMeshToSubmeshBary(o.m0, o.m0LvBridge, "tv", "tm");
o.m0RvBridge = mapCoordinateFromMeshToSubmeshBary(o.m0, o.m0RvBridge, "tv", "tm");

% visualize mapped coordinates
if o.cfg.exportLevel > 1
    vtkWrite(o.m0Ventricles.vol, sprintf('%sventricularmeshTmTv.vtu', o.cfg.outPrefix));
    vtkWrite(o.m0LvBridge.vol, sprintf('%slvBridgeTmTv.vtu', o.cfg.outPrefix));
    vtkWrite(o.m0RvBridge.vol, sprintf('%srvBridgeTmTv.vtu', o.cfg.outPrefix));
end

% map given boundary surfaces to surface of volume mesh and annotate classes
[o.m0LvBridge] = mapSurfaceClassesOnToBridge(o.m0LvBridge, lvOTBoundary,lvITBoundary, o.lvBridgeBoundaryFw, o.lvBridgeBoundarySeptum, o.cfg.outPrefix, o.cfg.mappingTol);
[o.m0RvBridge] = mapSurfaceClassesOnToBridge(o.m0RvBridge, rvOTBoundary,rvITBoundary, o.rvBridgeBoundaryFw, o.rvBridgeBoundarySeptum, o.cfg.outPrefix, o.cfg.mappingTol);

% Calculate gradient
P1 = double(o.m0RvBridge.vol.points);
C1 = double(o.m0RvBridge.vol.cells);
o.m0RvBridge.L = cotmatrix(P1, C1);
o.m0RvBridge.G = grad(P1, C1);

% calculate laplacian for FSW direction
idsBoundarySeptRv = o.m0RvBridge.surToVol(o.m0RvBridge.sur.pointData.class == 3);
idsBoundaryFwRv = o.m0RvBridge.surToVol(o.m0RvBridge.sur.pointData.class == 4);

ids = [idsBoundarySeptRv; idsBoundaryFwRv];
val = [zeros(size(idsBoundarySeptRv,1),1); ones(size(idsBoundaryFwRv))];
abLaplaceBridgeRv = solveLaplace(o.m0RvBridge.L, ids, val, o.cfg.tol, o.cfg.maxit);

o.m0RvBridge.vol.pointData.abLaplaceBridgeRv = single(abLaplaceBridgeRv);
vtkWrite(o.m0RvBridge.vol, sprintf('%sabLaplaceBridgeRvAbLaplacaceRt.vtu', o.cfg.outPrefix));

% calculate laplacian for anterior-posterior direction
% define boundary ids
idsBoundaryTvRv = o.m0RvBridge.surToVol(o.m0RvBridge.sur.pointData.class == 1);
idsBoundaryPvRv = o.m0RvBridge.surToVol(o.m0RvBridge.sur.pointData.class == 2);

ids = [idsBoundaryTvRv; idsBoundaryPvRv];
val = [zeros(size(idsBoundaryTvRv,1),1); ones(size(idsBoundaryPvRv))];
o.rtLaplaceBridgeRv = solveLaplace(o.m0RvBridge.L, ids, val, o.cfg.tol, o.cfg.maxit);
o.m0RvBridge.vol.pointData.rtLaplaceBridgeRv = o.rtLaplaceBridgeRv;

% set tolerance for normalized gradient field
tolGradField = 1e-8;

% calculating rotational coordinate
tmGrad = normalizedGradField(o.m0RvBridge.G, o.m0RvBridge.vol.pointData.tm, tolGradField, true, o.m0RvBridge.vol.points, o.m0RvBridge.vol.cells);
abLaplaceGrad = normalizedGradField(o.m0RvBridge.G, abLaplaceBridgeRv, tolGradField, true,  o.m0RvBridge.vol.points, o.m0RvBridge.vol.cells);
rtGrad = cross(tmGrad, abLaplaceGrad);
o.GrvBridge = grad(double(o.m0RvBridge.vol.points), double(o.m0RvBridge.vol.cells));

% Define anterior and posterior boundaries
idsPostrtrv = idsBoundaryTvRv;
idsAnttrtrv = idsBoundaryPvRv;

% Define tolerance for calculating TrajectProjection
tolTrajectProj = 1e-8;

drvbridgePost = solveTrajectDist(o.GrvBridge, rtGrad, idsPostrtrv, zeros(size(idsPostrtrv)), tolTrajectProj, o.cfg.maxit);
drvbridgeAnt = solveTrajectDist(o.GrvBridge, -rtGrad, idsAnttrtrv, zeros(size(idsAnttrtrv)), tolTrajectProj, o.cfg.maxit);
rtTrajectDistSept = drvbridgePost./(drvbridgeAnt+drvbridgePost);

% Scale PA from 1 to 1.5
rtSept = 1 + 0.5 * (1-rtTrajectDistSept);

o.m0RvBridge.vol.pointData.rtSin = sin(2*pi*rtSept);
o.m0RvBridge.vol.pointData.rtCos = cos(2*pi*rtSept);

o.m0RvBridge.vol.cellData.rtGrad = rtGrad;

if o.cfg.exportLevel > 1

    o.m0RvBridge.vol.pointData.drvbridgePost = single(drvbridgePost);
    vtkWrite(o.m0RvBridge.vol, sprintf('%srvDrvbridgePost.vtu', o.cfg.outPrefix));

    o.m0RvBridge.vol.pointData.drvbridgeAnt = single(drvbridgeAnt);
    vtkWrite(o.m0RvBridge.vol, sprintf('%srvDrvbridgeAnt.vtu', o.cfg.outPrefix));
end

o.m0RvBridge.vol.pointData.rt = rtSept;


% Compute anterior posterior in LV intervalvular regions
P1 = double(o.m0LvBridge.vol.points);
C1 = double(o.m0LvBridge.vol.cells);
o.m0LvBridge.L = cotmatrix(P1, C1);
o.m0LvBridge.G = grad(P1, C1);

% calculate laplacian for ab direction
idsBoundarySeptLv = o.m0LvBridge.surToVol(o.m0LvBridge.sur.pointData.class == 3);
idsBoundaryFwLv = o.m0LvBridge.surToVol(o.m0LvBridge.sur.pointData.class == 4);

ids = [idsBoundarySeptLv; idsBoundaryFwLv];
val = [zeros(size(idsBoundarySeptLv,1),1); ones(size(idsBoundaryFwLv))];
abLaplaceBridgeLv = solveLaplace(o.m0LvBridge.L, ids, val, o.cfg.tol, o.cfg.maxit);

o.m0LvBridge.vol.pointData.abLaplaceBridgeLv= single(abLaplaceBridgeLv);
vtkWrite(o.m0LvBridge.vol, sprintf('%sabLaplaceBridgeLv.vtu', o.cfg.outPrefix));

% calculate laplacian for pa direction
idsBoundaryMvLv = o.m0LvBridge.surToVol(o.m0LvBridge.sur.pointData.class == 1);
idsBoundaryAvLv = o.m0LvBridge.surToVol(o.m0LvBridge.sur.pointData.class == 2);

ids = [idsBoundaryMvLv; idsBoundaryAvLv];
val = [zeros(size(idsBoundaryMvLv,1),1); ones(size(idsBoundaryAvLv))];

o.rtLaplaceBridgeLv= solveLaplace(o.m0LvBridge.L, ids, val, o.cfg.tol, o.cfg.maxit);
o.m0LvBridge.vol.pointData.rtLaplaceBridgeLv= o.rtLaplaceBridgeLv;

tmGrad = normalizedGradField(o.m0LvBridge.G , o.m0LvBridge.vol.pointData.tm, tolGradField, true, o.m0LvBridge.vol.points, o.m0LvBridge.vol.cells);
abLaplaceGrad = normalizedGradField(o.m0LvBridge.G, abLaplaceBridgeLv, tolGradField, true,  o.m0LvBridge.vol.points, o.m0LvBridge.vol.cells);

rtGrad = cross(tmGrad, abLaplaceGrad);

o.GlvBridge = grad(double(o.m0LvBridge.vol.points), double(o.m0LvBridge.vol.cells));

idsPostrtlv = idsBoundaryMvLv;
idsAnttrtlv = idsBoundaryAvLv;

dlvbridgePost = solveTrajectDist(o.GlvBridge, rtGrad, idsPostrtlv, zeros(size(idsPostrtlv)), tolTrajectProj, o.cfg.maxit);
dlvbridgeAnt = solveTrajectDist(o.GlvBridge, -rtGrad, idsAnttrtlv, zeros(size(idsAnttrtlv)), tolTrajectProj, o.cfg.maxit);
rtTrajectDistSept = dlvbridgePost./(dlvbridgeAnt+dlvbridgePost);

% Scale PA > 1
rtSept = 1 + 0.5 * (1-rtTrajectDistSept);

o.m0LvBridge.vol.pointData.rtSin = sin(2*pi*rtSept);
o.m0LvBridge.vol.pointData.rtCos = cos(2*pi*rtSept);

o.m0LvBridge.vol.cellData.rtGrad = rtGrad;

if o.cfg.exportLevel > 1

    o.m0LvBridge.vol.pointData.dlvbridgePost = single(dlvbridgePost);
    vtkWrite(o.m0LvBridge.vol, sprintf('%slvDlvbridgePost.vtu', o.cfg.outPrefix));

    o.m0LvBridge.vol.pointData.dlvbridgeAnt = single(dlvbridgeAnt);
    vtkWrite(o.m0LvBridge.vol, sprintf('%slvDlvbridgeAnt.vtu', o.cfg.outPrefix));

end

%issue
o.m0LvBridge.vol.pointData.rt = rtSept;

o.printStatus(sprintf('%.1f seconds\n', toc-t), true);
o.available.rotationalbridges = true;

end
