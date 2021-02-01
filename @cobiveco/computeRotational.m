function computeRotational(o)

if ~o.available.heartAxesAndApex
    o.computeHeartAxesAndApex;
end
if ~o.available.mesh2
    o.prepareMesh2;
end
o.printStatus('Computing rotational coordinate...');
t = toc;

ucl = unique(o.m2.vol.cellData.class);
idsRidge = intersect(o.m2.vol.cells(o.m2.vol.cellData.class==ucl(1),:), o.m2.vol.cells(o.m2.vol.cellData.class==ucl(2),:));

septCurve = vtkDeleteDataArrays(vtkAppendPolyData({o.septCurveAnt, o.septCurvePost}));
mappedIds = knnsearch(septCurve.points, o.m2.vol.points(idsRidge,:), 'NSMethod','kdtree');
idsAnt = idsRidge(mappedIds <= size(o.septCurveAnt.points,1));
idsPost = setdiff(idsRidge, idsAnt);

ridge = repmat(uint8(0),size(o.m2.vol.points,1),1);
ridge(idsAnt) = 1;
ridge(idsPost) = 2;

apexCells = o.m2.vol.cells(any(ridge(o.m2.vol.cells)==1,2) & any(ridge(o.m2.vol.cells)==2,2),:);
idsApex = unique(apexCells(ridge(apexCells)==1));
apex = repmat(uint8(0),size(o.m2.vol.points,1),1);
apex(idsApex) = 1;

idsBase = o.m2.surToVol(o.m2.sur.pointData.class==1);

ids = [idsApex; idsBase];
val = [zeros(size(idsApex)); ones(size(idsBase))];
o.m2.abLaplace = solveLaplace(o.m2.L, ids, val, o.cfg.tol, o.cfg.maxit);

cellIdsSept = find(o.m2.vol.cellData.class==ucl(1));
cellIdsFree = find(o.m2.vol.cellData.class==ucl(2));

tmp = o.m2.vol;
tmp.pointData.ids = int32(1:size(tmp.points,1))';
sept = vtkThreshold(tmp, 'cells', 'class', double([ucl(1) ucl(1)]));
free = vtkThreshold(tmp, 'cells', 'class', double([ucl(2) ucl(2)]));

sept.pointData.ridge = ridge(sept.pointData.ids);
free.pointData.ridge = ridge(free.pointData.ids);

% make sure the epicardial part of the septal region has only non-zero ridge values
% by overwriting zero ridge values by their nearest-neighbor non-zero ridge values
idsSeptEpi = find(ismember(sept.pointData.ids, o.m2.surToVol(o.m2.sur.pointData.class==2)));
idsSeptTarget = idsSeptEpi(sept.pointData.ridge(idsSeptEpi)==0);
idsSeptSource = setdiff(idsSeptEpi, idsSeptTarget);
idsSeptSource = idsSeptSource(knnsearch(sept.points(idsSeptSource,:), sept.points(idsSeptTarget,:)));
sept.pointData.ridge(idsSeptTarget) = sept.pointData.ridge(idsSeptSource);
ridge(sept.pointData.ids(idsSeptTarget)) = sept.pointData.ridge(idsSeptSource);

idsSeptAnt = find(sept.pointData.ridge==1);
idsSeptPost = find(sept.pointData.ridge==2);
idsFreeAnt = find(free.pointData.ridge==1);
idsFreePost = find(free.pointData.ridge==2);

ucl = unique(o.m1.vol.cellData.class);
idsRight = o.m1.vol.cells(o.m1.vol.cellData.class==ucl(1),:);
tmFlipped = o.m1.tm;
tmFlipped(idsRight) = -tmFlipped(idsRight);
M12 = baryInterpMat(o.m1.vol.points, o.m1.vol.cells, o.m2.vol.points);
tmFlipped = M12 * tmFlipped;

tmGrad = normalizedGradField(o.m2.G, tmFlipped, o.cfg.tol, true, o.m2.vol.points, o.m2.vol.cells);
abLaplaceGrad = normalizedGradField(o.m2.G, o.m2.abLaplace, o.cfg.tol, true, o.m2.vol.points, o.m2.vol.cells);

rtGrad = cross(tmGrad, abLaplaceGrad);

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

rtTrajectDist = NaN(size(o.m2.vol.points,1),1);
rtTrajectDist(free.pointData.ids) = rtTrajectDistFree;
rtTrajectDist(sept.pointData.ids) = rtTrajectDistSept;

o.m2.rt = NaN(size(o.m2.vol.points,1),1);
o.m2.rt(free.pointData.ids) = rtFree;
o.m2.rt(sept.pointData.ids) = rtSept;
o.m2.rt = min(max(o.m2.rt,0),1);

o.m2.rtSin = sin(2*pi*o.m2.rt);
o.m2.rtCos = cos(2*pi*o.m2.rt);

o.m0.rtSin = min(max(o.m2.M*o.m2.rtSin,-1),1);
o.m0.rtCos = min(max(o.m2.M*o.m2.rtCos,-1),1);
o.m0.rt = atan2(o.m0.rtSin,o.m0.rtCos)/(2*pi);
o.m0.rt(o.m0.rt<0) = o.m0.rt(o.m0.rt<0)+1;
o.result.pointData.rtSin = single(o.m0.rtSin);
o.result.pointData.rtCos = single(o.m0.rtCos);
o.result.pointData.rt = single(o.m0.rt);

if o.cfg.exportLevel > 1
    o.m2.debug.pointData.ridge = ridge;
    o.m2.debug.pointData.apex = apex;
    o.m2.debug.pointData.abLaplace = single(o.m2.abLaplace);
    o.m2.debug.pointData.rtTrajectDist = single(rtTrajectDist);
    o.m2.debug.pointData.rt = single(o.m2.rt);
    o.m2.debug.pointData.rtSin = single(o.m2.rtSin);
    o.m2.debug.pointData.rtCos = single(o.m2.rtCos);
    vtkWrite(o.m2.debug, sprintf('%sdebug2.vtu', o.cfg.outPrefix));
    
    o.m0.debug.pointData.rt = single(o.m0.rt);
    o.m0.debug.pointData.rtSin = single(o.m0.rtSin);
    o.m0.debug.pointData.rtCos = single(o.m0.rtCos);
    vtkWrite(o.m0.debug, sprintf('%sdebug0.vtu', o.cfg.outPrefix));
end

o.printStatus(sprintf('%.1f seconds\n', toc-t), true);
o.available.rotational = true;

end