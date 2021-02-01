function computeHeartAxesAndApex(o)

if ~o.available.mesh1
    o.prepareMesh1;
end
o.printStatus('Computing heart axes and apex point...');
t = toc;

% point ids (wrt o.m1.vol) of septal surface and septal curve on epicaridum
idsEpi = o.m1.surToVol(o.m1.sur.pointData.class==2);
ucl = unique(o.m1.vol.cellData.class);
idsSeptSur = intersect(o.m1.vol.cells(o.m1.vol.cellData.class==ucl(1),:), o.m1.vol.cells(o.m1.vol.cellData.class==ucl(2),:));
idsSeptCurve = intersect(idsEpi, idsSeptSur);

% septCurve: line segments of septal curve, extracted from o.m1.sur
[~,idsSeptCurveSur] = ismember(idsSeptCurve, o.m1.surToVol);
pSeptCurve = double(o.m1.sur.points(idsSeptCurveSur,:));
cSeptCurve = o.m1.sur.cells';
ind = ismember(cSeptCurve, idsSeptCurveSur);
ind(:,sum(ind,1)~=2) = 0;
cSeptCurve = unique(reshape(cSeptCurve(ind),2,[])', 'rows');
cSeptCurve = changem(cSeptCurve, 1:numel(idsSeptCurveSur), idsSeptCurveSur);
septCurve = vtkCreateStruct(pSeptCurve, cSeptCurve);

% lvCenter: center of LV endocardium 
endoLv = vtkDataSetSurfaceFilter(vtkThreshold(o.m1.sur, 'points', 'class', [3 3]));
lvCenter = vtkCenterOfArea(endoLv);

% rvCenter: center of RV endocardium 
endoRv = vtkDataSetSurfaceFilter(vtkThreshold(o.m1.sur, 'points', 'class', [4 4]));
rvCenter = vtkCenterOfArea(endoRv);

% point coords of base
pBase = double(o.m1.sur.points(o.m1.sur.pointData.class==1,:));

% longAx: vector that is on average most orthogonal to the surface normals of the LV endocardium
longAx = computeLongAxis(endoLv, lvCenter-mean(pBase,1));

% compute triangle centroids and areas of septSur
pSeptSur = double(o.m1.vol.points(idsSeptSur,:));
cSeptSur = o.m1.vol.cells';
ind = ismember(cSeptSur, idsSeptSur');
ind(:,sum(ind,1)~=3) = 0;
cSeptSur = unique(reshape(cSeptSur(ind),3,[])', 'rows');
cSeptSur = changem(cSeptSur, 1:numel(idsSeptSur), idsSeptSur);
centroidsSeptSur = squeeze(mean(reshape(pSeptSur(cSeptSur,:),[],size(cSeptSur,2),size(pSeptSur,2)),2));
areasSeptSur = doublearea(pSeptSur, cSeptSur)/2;

% septVec: vector most orthogonal to the septal surface,
% determined as the last principle component of septal surface centroids
pc = pca(centroidsSeptSur, 'Weights', areasSeptSur./mean(areasSeptSur));
septVec = pc(:,end)';
% make sure septVec is directed from left to right
d = (rvCenter - lvCenter) * septVec';
septVec = sign(d) * septVec;

% antPostVec: vector pointing from anterior to posterior
antPostVec = cross(longAx, septVec);
d = (centroidsSeptSur - repmat(lvCenter, size(centroidsSeptSur,1), 1)) * antPostVec';
maskSeptSur = d > prctile(d,o.cfg.truncSeptSur(1)) & d < prctile(d,100-o.cfg.truncSeptSur(2));

% truncSeptCenter: center of truncated septal surface
centroidsTruncSept = centroidsSeptSur(maskSeptSur,:);
areasTruncSept = areasSeptSur(maskSeptSur);
truncSeptCenter = sum(repmat(areasTruncSept,1,3) .* centroidsTruncSept, 1) / sum(areasTruncSept);

% leftRightAx: vector pointing from left to right,
% determined as the last principle component of truncated septal surface centroids
pc = pca(centroidsTruncSept, 'Weights', areasTruncSept./mean(areasTruncSept));
leftRightAx = pc(:,end)';
% make sure leftRightAx is directed from left to right
d = (truncSeptCenter - lvCenter) * leftRightAx';
leftRightAx = sign(d) * leftRightAx;
% orthogonalize leftRightAx wrt longAx
leftRightAx = leftRightAx - leftRightAx*(longAx'*longAx);

% antPostAx: vector pointing from anterior to posterior,
antPostAx = cross(longAx, leftRightAx);

% center: global center point,
% determined by projecting lvCenter onto the truncated septal surface
d = (truncSeptCenter - lvCenter) * leftRightAx';
center = lvCenter + d .* leftRightAx;

% pApexCurve: points of septal curve below the center point
d = (pSeptCurve - repmat(center, size(pSeptCurve,1), 1)) * longAx';
ind1 = find(d>0);
pApexCurve = pSeptCurve(ind1,:);

% apex point: point in pApexCurve that is closest to the line parallel to
% longAx and passing through center
d = (repmat(center, size(pApexCurve,1), 1) - pApexCurve) * longAx';
pApexCurveProj = pApexCurve + d .* repmat(longAx, size(pApexCurve,1), 1);
dist = sqrt(sum((pApexCurveProj - center).^2, 2));
[~,ind2] = min(dist);

% split septCurve at apex point into anterior and posterior part
ind3 = ind1(ind2);
septCurveSplit = septCurve;
ind4 = any(ismember(septCurveSplit.cells, ind3),2);
septCurveSplit.cells(ind4,:) = [];
septCurveSplit.cellTypes(ind4,:) = [];
septCurveSplit = vtkConnectivityFilter(septCurveSplit);
isovals = double(unique(septCurveSplit.cellData.RegionId));
numCellsPerRegion = NaN(size(isovals));
for i = 1:numel(isovals)
    numCellsPerRegion(i) = numel(find(septCurveSplit.cellData.RegionId==isovals(i)));
end
[~,ind5] = sort(numCellsPerRegion, 'descend');
isovals = isovals(ind5(1:2));
o.septCurveAnt = vtkDeleteDataArrays(vtkThreshold(septCurveSplit, 'cells', 'RegionId', [isovals(1) isovals(1)]));
o.septCurvePost = vtkDeleteDataArrays(vtkThreshold(septCurveSplit, 'cells', 'RegionId', [isovals(2) isovals(2)]));
% make sure o.septCurveAnt and o.septCurvePost are not interchanged
if (mean(o.septCurvePost.points,1)-center)*antPostAx' < 0
    tmp = o.septCurveAnt;
    o.septCurveAnt = o.septCurvePost;
    o.septCurvePost = tmp;
end

% R: rotation matrix to align heart axes with x,y,z
xvec = -leftRightAx';
zvec = -longAx';
R1 = makehgtform('translate', -center);
R2 = makehgtform('axisrotate', cross([0;0;1], zvec), -acos([0;0;1]'*zvec));
R3 = makehgtform('zrotate', -acos([1;0;0]'*(R2(1:3,1:3)*xvec)));
o.R = R3*R2*R1;

if o.cfg.exportLevel > 2
    vtkWrite(o.septCurveAnt, sprintf('%sseptCurveAnt.vtp', o.cfg.outPrefix));
    vtkWrite(o.septCurvePost, sprintf('%sseptCurvePost.vtp', o.cfg.outPrefix));
    
    cSeptSurTrunc = cSeptSur(maskSeptSur,:);
    pSeptSurTrunc = pSeptSur(unique(cSeptSurTrunc),:);
    cSeptSurTrunc = changem(cSeptSurTrunc, 1:size(pSeptSurTrunc,1), unique(cSeptSurTrunc));
    vtkWrite(vtkCreateStruct(pSeptSurTrunc, cSeptSurTrunc), sprintf('%sseptSurfaceTrunc.vtp', o.cfg.outPrefix));
    
    d = norm(center-lvCenter);
    axes.points = [center; center-d*leftRightAx; center+2.5*d*longAx; center+2*d*antPostAx];
    axes.cells = int32([1 2; 1 3; 1 4]);
    axes.cellTypes = uint8([3; 3; 3]);
    axes.cellData.axis = uint8([1; 2; 3]);
    vtkWrite(axes, sprintf('%saxes.vtp', o.cfg.outPrefix));
end

o.printStatus(sprintf('%.1f seconds\n', toc-t), true);
o.available.heartAxesAndApex = true;

end