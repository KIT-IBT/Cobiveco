function computeApicobasal(o)

if ~o.available.splitMesh
    o.splitMesh;
end
o.printStatus('Computing apicobasal coordinate...');
t = toc;

idsBase = o.bv.surToVol(o.bv.sur.pointData.class==1);
idEpiApex = o.bv.surToVol(o.bv.sur.pointData.class==5);

idsEpiRv = intersect(o.bv.surToVol(o.bv.sur.pointData.class==2), find(o.bv.tv==1));
[~,ind] = min(sum((o.bv.vol.points(idsEpiRv,:)-repmat(o.bv.vol.points(idEpiApex,:),numel(idsEpiRv),1)).^2,2));
idEpiApex = idsEpiRv(ind);
o.bv.sur.pointData.class(o.bv.sur.pointData.class==5) = 2;
o.bv.sur.pointData.class(o.bv.surToVol==idEpiApex) = 5;

idsEndoLv = o.bv.surToVol(o.bv.sur.pointData.class==3);
[~,ind] = min(sum((o.bv.vol.points(idsEndoLv,:)-repmat(o.bv.vol.points(idEpiApex,:),numel(idsEndoLv),1)).^2,2));
idEndoLvApex = idsEndoLv(ind);

idsEndoRv = o.bv.surToVol(o.bv.sur.pointData.class==4);
[~,ind] = min(sum((o.bv.vol.points(idsEndoRv,:)-repmat(o.bv.vol.points(idEpiApex,:),numel(idsEndoRv),1)).^2,2));
idEndoRvApex = idsEndoRv(ind);

o.lvApexVec = o.bv.vol.points(idEndoLvApex,:)-o.bv.vol.points(idEpiApex,:);
o.rvApexVec = o.bv.vol.points(idEndoRvApex,:)-o.bv.vol.points(idEpiApex,:);

d1 = sqrt(sum(cross(o.bv.vol.points-o.bv.vol.points(idEpiApex,:), repmat(o.lvApexVec,size(o.bv.vol.points,1),1)).^2,2)/sum(o.lvApexVec.^2,2));
d2 = (o.bv.vol.points-o.bv.vol.points(idEpiApex,:))*o.lvApexVec'/norm(o.lvApexVec);
d3 = o.lvApexVec*o.lvApexVec'/norm(o.lvApexVec);
idsApex = find(d1<1.5*o.meanEdgLen & d2>-o.meanEdgLen & d2<d3+o.meanEdgLen);

d1 = sqrt(sum(cross(o.bv.vol.points-o.bv.vol.points(idEpiApex,:), repmat(o.rvApexVec,size(o.bv.vol.points,1),1)).^2,2)/sum(o.rvApexVec.^2,2));
d2 = (o.bv.vol.points-o.bv.vol.points(idEpiApex,:))*o.rvApexVec'/norm(o.rvApexVec);
d3 = o.rvApexVec*o.rvApexVec'/norm(o.rvApexVec);
idsApex = unique([idsApex; find(d1<1.5*o.meanEdgLen & d2>-o.meanEdgLen & d2<d3+o.meanEdgLen)]);

ids = [idsApex; idsBase];
val = [zeros(size(idsApex)); ones(size(idsBase))];
o.bv.abLaplace = solveLaplace(o.bv.L, ids, val, o.cfg.tol, o.cfg.maxit);

idsStart = find(o.bv.sur.pointData.class==1);
idEnd = find(o.bv.sur.pointData.class==5);
D = perform_fast_marching_mesh(double(o.bv.sur.points), double(o.bv.sur.cells), idsStart);
geodesicPoints = compute_geodesic_mesh(D, double(o.bv.sur.points), double(o.bv.sur.cells), idEnd, struct('method','continuous'))';
tr = vtkToTriangulation(o.bv.sur);
geodesic = vtkCreateStruct(geodesicPoints);
geodesic.pointData.normals = tr.vertexNormal(double(vtkMapPointIds(o.bv.sur, geodesic)));
M = vtkBarycentricInterpMat(o.bv.sur, geodesic, o.meanEdgLen);
geodesicAb = M * o.bv.abLaplace(vtkMapPointIds(o.bv.vol, o.bv.sur));
nanInd = isnan(geodesicAb);
geodesicAb(nanInd) = [];
geodesic.points(nanInd,:) = [];
geodesic.pointData.normals(nanInd,:) = [];
geodesic.pointData.abLaplace = geodesicAb;
geodesicDist = [0; cumsum(sqrt(sum(diff(geodesicPoints,1,1).^2,2)))];
geodesicDist = geodesicDist/max(geodesicDist);
geodesicDist(nanInd) = [];
[geodesicAb,ind] = unique(geodesicAb);
geodesicDist = geodesicDist(ind);
geodesicDist = [-1e-10; geodesicDist; 1+1e-10];
geodesicAb = [-1e-10; geodesicAb; 1+1e-10];
o.bv.ab = interp1(geodesicAb, geodesicDist, o.bv.abLaplace);

o.result.pointData.ab = single(o.bv.ab);

if o.cfg.exportLevel > 1
    vtkWrite(geodesic, sprintf('%sgeodesic.vtp', o.cfg.outPrefix));
end
if o.cfg.exportLevel > 2
    vtkWrite(o.result, sprintf('%sresult.vtu', o.cfg.outPrefix));
end

o.printStatus(sprintf('%.1f seconds\n', toc-t), true);
o.available.apicobasal = true;

end
