function computeRotationalVentricles(o)

% Computes the rotational coordinates of the object in the ventricles,
% i.e. the trajectory distance from posterior to anterior and back to posterior.
% The coordinate values are added in the fields "rt", "rtSin" and "rtCos" to the output objects.
%
% computeRotationalVentricles(o)
%
% Input:
%   o: instance of the class cobiveco
%
% Output:
%   o.m0, object of class cobiveco [struct] (For details see cobiveco class documentation)



if ~o.available.apicobasalbridges
    o.computeApicobasalBridges;
end

o.printStatus('Computing RT on ventricles.');
t = toc;

% map surface classes onto ventricular mesh
o.m0Ventricles.sur = vtkDataSetSurfaceFilter(vtkDeleteDataArrays(o.m0Ventricles.vol));
o.m0Ventricles.sur = vtkArrayMapperNearestNeighbor(o.m0.sur, o.m0Ventricles.sur);
o.m0Ventricles.surToVol = vtkMapPointIds(o.m0Ventricles.vol, o.m0Ventricles.sur);
o.m0Ventricles.meanEdgLen = mean(vtkEdgeLengths(o.m0Ventricles.vol));

P2 = double(o.m0Ventricles.vol.points);
C2 = double(o.m0Ventricles.vol.cells);
o.m0Ventricles.L = cotmatrix(P2, C2);
o.m0Ventricles.G = grad(P2, C2);

% Transfer classes - can be condensed as assigned twice
ventricularmeshParams = vtkArrayMapperNearestNeighbor(o.m1.vol, o.m0Ventricles.vol);
vtkWrite(o.m0Ventricles.vol, sprintf('%sventricularmeshclasses.vtu', o.cfg.outPrefix));

epiIds = o.m0Ventricles.surToVol(o.m0Ventricles.sur.pointData.class == 2);
ucl = unique(ventricularmeshParams.cellData.class);
septIds = intersect(ventricularmeshParams.cells(ventricularmeshParams.cellData.class == ucl(1),:), ventricularmeshParams.cells(ventricularmeshParams.cellData.class == ucl(2),:));
epiIds = setdiff(epiIds, septIds);
ids = [epiIds; septIds];
val = [zeros(size(epiIds)); ones(size(septIds))];
o.m0Ventricles.ridgeLaplace = solveLaplace(o.m0Ventricles.L, ids, val, o.cfg.tol, o.cfg.maxit);

o.m0Ventricles.vol.pointData.ridgeLaplace = o.m0Ventricles.ridgeLaplace;
vtkWrite(ventricularmeshParams, sprintf('%sventricularmeshParams.vtu', o.cfg.outPrefix));
o.m0Ventricles.vol.pointData.ridgeLaplace = o.m0Ventricles.ridgeLaplace;
vtkWrite(o.m0Ventricles.vol, sprintf('%sventricularmeshPrep.vtu', o.cfg.outPrefix));


% define Epi Ids
idsEpi = zeros(size(o.m0Ventricles.vol.points,1),1);
idsEpi(epiIds) = 1;
o.m0Ventricles.vol.pointData.idsEpi = idsEpi;
idsSept = zeros(size(o.m0Ventricles.vol.points,1),1);
idsSept(septIds) = 1;
o.m0Ventricles.vol.pointData.idsSept = idsSept;
vtkWrite(o.m0Ventricles.vol, sprintf('%sventricularmeshPrep.vtu', o.cfg.outPrefix));

% calculate tv laplace and remesh
idsEndoLv = o.m0Ventricles.surToVol(o.m0Ventricles.sur.pointData.class == 3);
idsEndoRv = o.m0Ventricles.surToVol(o.m0Ventricles.sur.pointData.class == 4);

ids = [idsEndoLv; idsEndoRv];
val = [zeros(size(idsEndoLv)); ones(size(idsEndoRv))];
o.m0Ventricles.tvLaplace = solveLaplace(o.m0Ventricles.L, ids, val, o.cfg.tol, o.cfg.maxit);
o.m0Ventricles.tv = round(o.m0Ventricles.tvLaplace);
o.m0Ventricles.vol.pointData.tv = single(o.m0Ventricles.tv);
o.m0Ventricles.vol.pointData.tvLaplace = single(o.m0Ventricles.tvLaplace);

% call remeshing function equals mesh1 for tv separation
cfg = o.cfg;
o.m1Ventricles = prepareSeptalSurface(o.m0Ventricles, cfg);
vtkWrite(o.m1Ventricles.vol, sprintf('%snewmesventricularmeshRemeshedp.vtu', o.cfg.outPrefix));

epiIds = o.m1Ventricles.surToVol(o.m1Ventricles.sur.pointData.class == 2);
ucl = unique(o.m1Ventricles.vol.cellData.class);
septIds = intersect(o.m1Ventricles.vol.cells(o.m1Ventricles.vol.cellData.class == ucl(1),:), o.m1Ventricles.vol.cells(o.m1Ventricles.vol.cellData.class == ucl(2),:));
idsEpiNoBase = setdiff(epiIds, septIds);

ids = [idsEpiNoBase; septIds];
val = [zeros(size(idsEpiNoBase)); ones(size(septIds))];
o.m1Ventricles.ridgeLaplace = solveLaplace(o.m1Ventricles.L, ids, val, o.cfg.tol, o.cfg.maxit);
o.m1Ventricles.vol.pointData.ridgeLaplace = o.m1Ventricles.ridgeLaplace;
idsEpi = zeros(size(o.m1Ventricles.vol.points,1),1);
idsEpi(idsEpiNoBase) = 1;
o.m1Ventricles.vol.pointData.idsEpi = idsEpi;
idsSept = zeros(size(o.m1Ventricles.vol.points,1),1);
idsSept(septIds) = 1;
o.m1Ventricles.vol.pointData.idsSept = idsSept;
vtkWrite(o.m1Ventricles.vol, sprintf('%snewmesventricularmeshRemeshedLaplacian.vtu', o.cfg.outPrefix));

original = o.m0;
cfg = o.cfg;
t = toc;

% remesh #2
o.m2Ventricles = prepareMeshVentricularVolume(o.m1Ventricles, original, cfg);

% Compute rotational
ucl = unique(o.m2Ventricles.vol.cellData.class);
idsRidge = intersect(o.m2Ventricles.vol.cells(o.m2Ventricles.vol.cellData.class == ucl(1),:), o.m2Ventricles.vol.cells(o.m2Ventricles.vol.cellData.class == ucl(2),:));

% Calculate septal curve on ventricular mesh
septCurve = vtkDeleteDataArrays(vtkAppendPolyData({o.septCurveAnt, o.septCurvePost}));
mappedIds = knnsearch(septCurve.points, o.m2Ventricles.vol.points(idsRidge,:), 'NSMethod','kdtree');
idsAnt = idsRidge(mappedIds <= size(o.septCurveAnt.points,1));
idsPost = setdiff(idsRidge, idsAnt);

ridge = repmat(uint8(0),size(o.m2Ventricles.vol.points,1),1);
ridge(idsAnt) = 1;
ridge(idsPost) = 2;

% find apex on ventricle
apexCells = o.m2Ventricles.vol.cells(any(ridge(o.m2Ventricles.vol.cells) == 1,2) & any(ridge(o.m2Ventricles.vol.cells) == 2,2),:);
idsApex = unique(apexCells(ridge(apexCells) == 1));
apex = repmat(uint8(0),size(o.m2Ventricles.vol.points,1),1);
apex(idsApex) = 1;

% Set as class instance! move to method
% set new boundary surfaces
surEpiBase = vtkRead(sprintf('%s_epi_base.ply', o.cfg.inPrefix));
surEpiNoBase = vtkRead(sprintf('%s_epi_no_base.ply', o.cfg.inPrefix));
surEndoLv = vtkRead(sprintf('%s_endo_lv.ply', o.cfg.inPrefix));
surEndoRv = vtkRead(sprintf('%s_endo_rv.ply', o.cfg.inPrefix));
o.m2Ventricles.meanEdgLen = mean(vtkEdgeLengths(o.m2Ventricles.vol));
boundaryLvBridgeSurAppended = vtkAppendPolyData({surEpiBase, surEpiNoBase, surEndoLv, surEndoRv, o.surMv,o.surTv, o.surAv, o.surPv,o.ventricularmeshRvBoundary,o.ventricularmeshLvBoundary});
cnp = cumsum([size(surEpiBase.points,1); size(surEpiNoBase.points,1);size(surEndoLv.points,1); size(surEndoRv.points,1);size(o.surMv.points,1); size(o.surTv.points,1); size(o.surAv.points,1); size(o.surPv.points,1); size(o.ventricularmeshRvBoundary.points,1); size(o.ventricularmeshLvBoundary.points,1)]);
ids = vtkMapPointIds(boundaryLvBridgeSurAppended, o.m2Ventricles.sur);
o.m2Ventricles.sur.pointData.newclass = repmat(uint8(0),size(o.m2Ventricles.sur.points,1),1);
o.m2Ventricles.sur.pointData.newclass(ids<=cnp(1)) = 1;               % EpiBase
o.m2Ventricles.sur.pointData.newclass(ids>cnp(1) & ids<=cnp(2)) = 2;  % surEpiNoBase
o.m2Ventricles.sur.pointData.newclass(ids>cnp(2) & ids<=cnp(3)) = 3;  % LV
o.m2Ventricles.sur.pointData.newclass(ids>cnp(3) & ids<=cnp(4)) = 4;  % Rv
o.m2Ventricles.sur.pointData.newclass(ids>cnp(4) & ids<=cnp(5)) = 5;  % MV
o.m2Ventricles.sur.pointData.newclass(ids>cnp(5) & ids<=cnp(6)) = 6;  % TV
o.m2Ventricles.sur.pointData.newclass(ids>cnp(6) & ids<=cnp(7)) = 7;  % AV
o.m2Ventricles.sur.pointData.newclass(ids>cnp(7) & ids<=cnp(8)) = 8;  % PV
o.m2Ventricles.sur.pointData.newclass(ids>cnp(8) & ids<=cnp(9)) = 9; % suro.m2VentriclesRv
o.m2Ventricles.sur.pointData.newclass(ids>cnp(9)) = 10; % suro.m2VentriclesLv
mappingDist = sqrt(sum((o.m2Ventricles.sur.points-boundaryLvBridgeSurAppended.points(ids,:)).^2,2));
o.m2Ventricles.sur.pointData.newclass(mappingDist > o.cfg.mappingTol*o.m2Ventricles.meanEdgLen) = 0;

o.m2Ventricles.vol.pointData.newsurClass = repmat(uint8(0),size(o.m2Ventricles.vol.points,1),1);
o.m2Ventricles.vol.pointData.newsurClass(o.m2Ventricles.surToVol) = o.m2Ventricles.sur.pointData.newclass;
vtkWrite(o.m2Ventricles.vol, sprintf('%so.m2Ventriclesclassesew.vtu', o.cfg.outPrefix));

idsBase = o.m2Ventricles.surToVol(o.m2Ventricles.sur.pointData.newclass == 5 |o.m2Ventricles.sur.pointData.newclass == 6 | o.m2Ventricles.sur.pointData.newclass == 7| o.m2Ventricles.sur.pointData.newclass == 8| o.m2Ventricles.sur.pointData.newclass == 9| o.m2Ventricles.sur.pointData.newclass == 10);
ids = [idsApex; idsBase];
val = [zeros(size(idsApex)); ones(size(idsBase))];
o.m2Ventricles.abLaplace = solveLaplace(o.m2Ventricles.L, ids, val, o.cfg.tol, o.cfg.maxit);
o.m2Ventricles.vol.pointData.abLaplace = o.m2Ventricles.abLaplace;
o.m2Ventricles.vol.pointData.ridge = ridge;
vtkWrite(o.m2Ventricles.vol, sprintf('%snewmeshablaplace.vtu', o.cfg.outPrefix));

cellIdsSept = find(o.m2Ventricles.vol.cellData.class == ucl(1));
cellIdsFree = find(o.m2Ventricles.vol.cellData.class == ucl(2));

tmp = o.m2Ventricles.vol;
tmp.pointData.ids = int32(1:size(tmp.points,1))';
sept = vtkThreshold(tmp, 'cells', 'class', double([ucl(1) ucl(1)]));
free = vtkThreshold(tmp, 'cells', 'class', double([ucl(2) ucl(2)]));

sept.pointData.ridge = ridge(sept.pointData.ids);
free.pointData.ridge = ridge(free.pointData.ids);

% make sure the epicardial part of the septal region has only non-zero ridge values
% by overwriting zero ridge values by their nearest-neighbor non-zero ridge values
idsSeptEpi = find(ismember(sept.pointData.ids, o.m2Ventricles.surToVol(o.m2Ventricles.sur.pointData.class == 2))); % add valve planes to boundary for epi!
idsSeptTarget = idsSeptEpi(sept.pointData.ridge(idsSeptEpi) == 0);
idsSeptSource = setdiff(idsSeptEpi, idsSeptTarget);
idsSeptSource = idsSeptSource(knnsearch(sept.points(idsSeptSource,:), sept.points(idsSeptTarget,:)));
sept.pointData.ridge(idsSeptTarget) = sept.pointData.ridge(idsSeptSource);
ridge(sept.pointData.ids(idsSeptTarget)) = sept.pointData.ridge(idsSeptSource);

idsSeptAnt = find(sept.pointData.ridge == 1);
idsSeptPost = find(sept.pointData.ridge == 2);
idsFreeAnt = find(free.pointData.ridge == 1);
idsFreePost = find(free.pointData.ridge == 2);

% mesh after first remeshing
ucl = unique(o.m1Ventricles.vol.cellData.class);
idsRight = o.m1Ventricles.vol.cells(o.m1Ventricles.vol.cellData.class == ucl(1),:);
% transfer tm
P1 = double(o.m1Ventricles.vol.points);
C1 = double(o.m1Ventricles.vol.cells);
o.m1Ventricles.L = cotmatrix(P1, C1);
o.m1Ventricles.G = grad(P1, C1);
o.m1Ventricles.massMat = massmatrix(P1, C1, 'voronoi');
Matrix = baryInterpMat(o.m1.vol.points,o.m1.vol.cells,o.m1Ventricles.vol.points);
o.m1Ventricles.tm = Matrix * o.m1.tm;
o.m1Ventricles.vol.pointData.tm = o.m1Ventricles.tm;

tmFlipped = o.m1Ventricles.tm;
tmFlipped(idsRight) = -tmFlipped(idsRight);
M12 = baryInterpMat(o.m1Ventricles.vol.points, o.m1Ventricles.vol.cells, o.m2Ventricles.vol.points);
tmFlipped = M12 * tmFlipped;

tmGrad = normalizedGradField(o.m2Ventricles.G, tmFlipped, o.cfg.tol, true, o.m2Ventricles.vol.points, o.m2Ventricles.vol.cells);
abLaplaceGrad = normalizedGradField(o.m2Ventricles.G, o.m2Ventricles.abLaplace, o.cfg.tol, true, o.m2Ventricles.vol.points, o.m2Ventricles.vol.cells);

rtGrad = cross(tmGrad, abLaplaceGrad);
o.m2Ventricles.vol.cellData.rtGrad = rtGrad;
GSept = grad(double(sept.points), double(sept.cells));
GFree = grad(double(free.points), double(free.cells));

dSeptPost = solveTrajectDist(GSept, rtGrad(cellIdsSept,:), idsSeptPost, zeros(size(idsSeptPost)), o.cfg.tol, o.cfg.maxit);
dSeptAnt = solveTrajectDist(GSept, -rtGrad(cellIdsSept,:), idsSeptAnt, zeros(size(idsSeptAnt)), o.cfg.tol, o.cfg.maxit);
rtTrajectDistSept = dSeptPost./(dSeptAnt+dSeptPost);
rtSept = 2/3 + 1/3 * (1-rtTrajectDistSept);

dFreePost = solveTrajectDist(GFree, rtGrad(cellIdsFree,:), idsFreePost, zeros(size(idsFreePost)), o.cfg.tol, o.cfg.maxit);
dFreeAnt = solveTrajectDist(GFree, -rtGrad(cellIdsFree,:), idsFreeAnt, zeros(size(idsFreeAnt)), o.cfg.tol, o.cfg.maxit);
rtTrajectDistFree = dFreePost./(dFreeAnt+dFreePost);
rtFree = 2/3 * rtTrajectDistFree;

rtTrajectDist = NaN(size(o.m2Ventricles.vol.points,1),1);
rtTrajectDist(free.pointData.ids) = rtTrajectDistFree;
rtTrajectDist(sept.pointData.ids) = rtTrajectDistSept;

o.m2Ventricles.rt = NaN(size(o.m2Ventricles.vol.points,1),1);
o.m2Ventricles.rt(free.pointData.ids) = rtFree;
o.m2Ventricles.rt(sept.pointData.ids) = rtSept;
o.m2Ventricles.rt = min(max(o.m2Ventricles.rt,0),1);

o.m2Ventricles.rtSin = sin(2*pi*o.m2Ventricles.rt);
o.m2Ventricles.rtCos = cos(2*pi*o.m2Ventricles.rt);

% Check rotational on ventricular volume
o.m2Ventricles.vol.pointData.rt = o.m2Ventricles.rt;
o.m2Ventricles.vol.pointData.rtSin = o.m2Ventricles.rtSin;
o.m2Ventricles.vol.pointData.rtCos = o.m2Ventricles.rtCos;
vtkWrite(o.m2Ventricles.vol, sprintf('%snewmeshvolvolRt.vtu', o.cfg.outPrefix));

% add to original mesh
o.m0.rtSin = min(max(o.m2Ventricles.M*o.m2Ventricles.rtSin,-1),1);
o.m0.rtCos = min(max(o.m2Ventricles.M*o.m2Ventricles.rtCos,-1),1);
o.m0.rt = atan2(o.m0.rtSin,o.m0.rtCos)/(2*pi);
o.m0.rt(o.m0.rt<0) = o.m0.rt(o.m0.rt<0)+1;
o.result.pointData.rtSin = single(o.m0.rtSin);
o.result.pointData.rtCos = single(o.m0.rtCos);
o.result.pointData.rt = single(o.m0.rt);
vtkWrite(o.result, sprintf('%sresultRt.vtu', o.cfg.outPrefix));


o.printStatus(sprintf('%.1f seconds\n', toc-t), true);
o.available.rotationalbridges = true;


end
