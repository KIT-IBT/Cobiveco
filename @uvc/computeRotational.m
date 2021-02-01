function computeRotational(o)

if ~o.available.apicobasal
    o.computeApicobasal;
end
o.printStatus('Computing rotational coordinate...');
t = toc;

%% boundaries of the LV-RV junction and first intermediate Laplace solution

junc = vtkFeatureEdges(vtkThreshold(o.lv.sur, 'cells', 'class', [0 0]), 1, 0, 0);
junc.pointData.abLaplace = o.bv.abLaplace(o.bv.surToVol(vtkMapPointIds(o.bv.sur, junc)));
juncEpi = vtkThreshold(junc, 'points', 'class', [2 2]);
juncEndo = vtkThreshold(junc, 'points', 'class', [6 6]);

juncEpi = vtkConnectivityFilter(juncEpi);
if numel(unique(juncEpi.pointData.RegionId)) == 1
    % default
    rid = double(mode(juncEpi.pointData.RegionId));
    juncEpi = vtkThreshold(juncEpi, 'points', 'RegionId', [rid rid]);
    abLaplaceSorted = sort(juncEpi.pointData.abLaplace);
    for k = 1:numel(abLaplaceSorted)
        juncEpi = vtkConnectivityFilter(vtkThreshold(juncEpi, 'points', 'abLaplace', [abLaplaceSorted(k)+o.cfg.tol inf]));
        [rid,~,subs] = unique(juncEpi.pointData.RegionId);
        if numel(rid) > 1
            freq = accumarray(subs,subs,[],@numel);
            [freq,sortInd] = sort(freq,'descend');
            if freq(2) > size(juncEpi.points,1)/4
                break;
            end
        end
    end
    if k == numel(abLaplaceSorted)
        error('juncEpi could not be splitted by thresholding abLaplace.');
    end
    rid1 = double(rid(sortInd(1)));
    rid2 = double(rid(sortInd(2)));
    juncEpiAnt = vtkThreshold(juncEpi, 'points', 'RegionId', [rid1 rid1]);
    juncEpiPost = vtkThreshold(juncEpi, 'points', 'RegionId', [rid2 rid2]);
else
    % fallback
    warning('Using fallback solution to split juncEpi.');
    [~,apexInd] = min(juncEpi.pointData.abLaplace);
    apexPoint = juncEpi.points(apexInd,:);
    basePoint = mean(o.lv.sur.points(o.lv.sur.pointData.class==1,:),1);
    leftRightVec = pca(juncEpi.points);
    leftRightVec = leftRightVec(:,end)';
    normalVec = cross(apexPoint-basePoint, leftRightVec);
    normalVec = normalVec/norm(normalVec);
    d = dot(juncEpi.points-apexPoint, repmat(normalVec,size(juncEpi.points,1),1), 2);
    juncEpi.pointData.d = d;
    juncEpiAnt = vtkThreshold(juncEpi, 'points', 'd', [-inf 0]);
    juncEpiPost = vtkThreshold(juncEpi, 'points', 'd', [0 inf]);
end

juncEndo = vtkConnectivityFilter(juncEndo);
if numel(unique(juncEndo.pointData.RegionId)) == 1
    % default
    rid = double(mode(juncEndo.pointData.RegionId));
    juncEndo = vtkThreshold(juncEndo, 'points', 'RegionId', [rid rid]);
    abLaplaceSorted = sort(juncEndo.pointData.abLaplace);
    for k = 1:numel(abLaplaceSorted)
        juncEndo = vtkConnectivityFilter(vtkThreshold(juncEndo, 'points', 'abLaplace', [abLaplaceSorted(k)+o.cfg.tol inf]));
        [rid,~,subs] = unique(juncEndo.pointData.RegionId);
        if numel(rid) > 1
            freq = accumarray(subs,subs,[],@numel);
            [freq,sortInd] = sort(freq,'descend');
            if freq(2) > size(juncEndo.points,1)/4
                break;
            end
        end
    end
    if k == numel(abLaplaceSorted)
        error('juncEndo could not be splitted by thresholding abLaplace.');
    end
    rid1 = double(rid(sortInd(1)));
    rid2 = double(rid(sortInd(2)));
    juncEndoAnt = vtkThreshold(juncEndo, 'points', 'RegionId', [rid1 rid1]);
    juncEndoPost = vtkThreshold(juncEndo, 'points', 'RegionId', [rid2 rid2]);
else
    % fallback
    warning('Using fallback solution to split juncEndo.');
    [~,apexInd] = min(juncEndo.pointData.abLaplace);
    apexPoint = juncEndo.points(apexInd,:);
    basePoint = mean(o.lv.sur.points(o.lv.sur.pointData.class==1,:),1);
    leftRightVec = pca(juncEndo.points);
    leftRightVec = leftRightVec(:,end)';
    normalVec = cross(apexPoint-basePoint, leftRightVec);
    normalVec = normalVec/norm(normalVec);
    d = dot(juncEndo.points-apexPoint, repmat(normalVec,size(juncEndo.points,1),1), 2);
    juncEndo.pointData.d = d;
    juncEndoAnt = vtkThreshold(juncEndo, 'points', 'd', [-inf 0]);
    juncEndoPost = vtkThreshold(juncEndo, 'points', 'd', [0 inf]);
end

baseToApex = mean(o.bv.sur.points(o.bv.sur.pointData.class==1,:),1)-o.bv.sur.points(o.bv.sur.pointData.class==5,:);
leftToRight = mean(o.bv.sur.points(o.bv.sur.pointData.class==3,:),1)-mean(o.bv.sur.points(o.bv.sur.pointData.class==4,:),1);
antToPost = cross(baseToApex, leftToRight);
antToPost = antToPost/norm(antToPost);

if (mean(juncEpiAnt.points,1)-mean(junc.points,1))*antToPost' > 0
    tmp = juncEpiAnt;
    juncEpiAnt = juncEpiPost;
    juncEpiPost = tmp;
end
if (mean(juncEndoAnt.points,1)-mean(junc.points,1))*antToPost' > 0
    tmp = juncEndoAnt;
    juncEndoAnt = juncEndoPost;
    juncEndoPost = tmp;
end

if o.cfg.exportLevel > 1
    vtkWrite(vtkDeleteDataArrays(juncEpiAnt), sprintf('%slvJuncEpiAnt.vtp', o.cfg.outPrefix));
    vtkWrite(vtkDeleteDataArrays(juncEpiPost), sprintf('%slvJuncEpiPost.vtp', o.cfg.outPrefix));
    vtkWrite(vtkDeleteDataArrays(juncEndoAnt), sprintf('%slvJuncEndoAnt.vtp', o.cfg.outPrefix));
    vtkWrite(vtkDeleteDataArrays(juncEndoPost), sprintf('%slvJuncEndoPost.vtp', o.cfg.outPrefix));
end

Mdl = KDTreeSearcher(o.lv.vol.points);
idsJuncEpiAnt = knnsearch(Mdl, juncEpiAnt.points);
idsJuncEpiPost = setdiff(knnsearch(Mdl, juncEpiPost.points), idsJuncEpiAnt);
idsJuncEndoAnt = setdiff(knnsearch(Mdl, juncEndoAnt.points), [idsJuncEpiAnt; idsJuncEpiPost]);
idsJuncEndoPost = setdiff(knnsearch(Mdl, juncEndoPost.points), [idsJuncEpiAnt; idsJuncEpiPost; idsJuncEndoAnt]);

ids = [idsJuncEpiAnt; idsJuncEndoAnt; idsJuncEpiPost; idsJuncEndoPost];
val = [ones(size([idsJuncEpiAnt; idsJuncEndoAnt])); -ones(size([idsJuncEpiPost; idsJuncEndoPost]))];
o.lv.lap1 = solveLaplace(o.lv.L, ids, val, o.cfg.tol, o.cfg.maxit);

%% zero and +-pi boundaries and second intermediate Laplace solution

zeroCross1.cells = o.lv.vol.cells(any(diff(o.lv.lap1(o.lv.vol.cells)<0,1,2)~=0,2),:);
zeroCross1.cellTypes = repmat(uint8(10), size(zeroCross1.cells,1), 1);
pointIds = unique(zeroCross1.cells(:));
zeroCross1.points = o.lv.vol.points(pointIds,:);
zeroCross1.pointData.ids = pointIds;
zeroCross1.pointData.abLaplace = o.bv.abLaplace(o.lv.ids(pointIds));
zeroCross1.pointData.lap1 = o.lv.lap1(pointIds);
zeroCross1.cells = changem(zeroCross1.cells, 1:size(zeroCross1.points,1), pointIds);
zeroCross1 = vtkConnectivityFilter(zeroCross1);
rid = double(mode(zeroCross1.pointData.RegionId));
zeroCross1 = vtkThreshold(zeroCross1, 'points', 'RegionId', [rid rid]);

% split at apical line into septal and free wall parts
zeroCross1Apex = vtkThreshold(zeroCross1, 'points', 'abLaplace', [-inf 0.4]);
zeroCross1Apex.pointData.ids = vtkMapPointIds(zeroCross1, zeroCross1Apex);
normalVec = pca(zeroCross1Apex.points);
normalVec = normalVec(:,end)';
apexVec = o.lvApexVec-o.lvApexVec*(normalVec'*normalVec);
posVec = o.bv.sur.points(o.bv.sur.pointData.class==5,:);
zeroCross1Apex.points = zeroCross1Apex.points-(zeroCross1Apex.points-posVec)*(normalVec'*normalVec);
d = sqrt(sum(cross(zeroCross1Apex.points-posVec, repmat(apexVec,size(zeroCross1Apex.points,1),1)).^2,2)/sum(apexVec.^2,2));
ids = zeroCross1Apex.pointData.ids(d<1.5*o.meanEdgLen);
zeroCross1.pointData.apex = zeros(size(zeroCross1.points,1),1);
zeroCross1.pointData.apex(ids) = 1;

zeroCross1 = vtkConnectivityFilter(vtkThreshold(zeroCross1, 'points', 'apex', [-inf 0]));
[rid,~,subs] = unique(zeroCross1.pointData.RegionId);
freq = accumarray(subs,subs,[],@numel);
[~,sortInd] = sort(freq,'descend');
rid1 = double(rid(sortInd(1)));
rid2 = double(rid(sortInd(2)));
zeroCross1Sept = vtkThreshold(zeroCross1, 'points', 'RegionId', [rid1 rid1]);
zeroCross1Free = vtkThreshold(zeroCross1, 'points', 'RegionId', [rid2 rid2]);

if (mean(zeroCross1Sept.points,1)-mean(zeroCross1.points,1))*leftToRight' > 0
    tmp = zeroCross1Sept;
    zeroCross1Sept = zeroCross1Free;
    zeroCross1Free = tmp;
end

idsZero = zeroCross1Sept.pointData.ids(zeroCross1Sept.pointData.lap1<0);
idsMinusPi = zeroCross1Free.pointData.ids(zeroCross1Free.pointData.lap1<0);
idsPlusPi = zeroCross1Free.pointData.ids(zeroCross1Free.pointData.lap1>0);

if o.cfg.exportLevel > 1
    vtkWrite(vtkDeleteDataArrays(vtkThreshold(vtkDataSetSurfaceFilter(zeroCross1Sept), 'points', 'lap1', [-inf 0])), sprintf('%slvZero.vtp', o.cfg.outPrefix));
    vtkWrite(vtkDeleteDataArrays(vtkThreshold(vtkDataSetSurfaceFilter(zeroCross1Free), 'points', 'lap1', [-inf 0])), sprintf('%slvMinusPi.vtp', o.cfg.outPrefix));
    vtkWrite(vtkDeleteDataArrays(vtkThreshold(vtkDataSetSurfaceFilter(zeroCross1Free), 'points', 'lap1', [0 inf])), sprintf('%slvPlusPi.vtp', o.cfg.outPrefix));
end

ids = [idsZero; idsMinusPi; idsPlusPi];
val = [zeros(size(idsZero)); repmat(-pi,size(idsMinusPi)); repmat(pi,size(idsPlusPi))];
o.lv.lap2 = solveLaplace(o.lv.L, ids, val, o.cfg.tol, o.cfg.maxit);

%% third and fourth intermediate Laplace solutions and +-pi/2 and +-pi/2.5 boundaries

meanEndoAnt  = mean(o.lv.lap2(vtkMapPointIds(o.lv.vol, juncEndoAnt)));
meanEndoPost = mean(o.lv.lap2(vtkMapPointIds(o.lv.vol, juncEndoPost)));
idsPlus = unique([ ...
    intersect(o.lv.surToVol(o.lv.sur.pointData.class==3), find(o.lv.lap2>meanEndoPost & o.lv.lap2<meanEndoAnt)); ... % endo
    o.lv.surToVol(o.lv.sur.pointData.class==6) ... % epi
    ]);
idsMinus = setdiff(o.lv.surToVol(o.lv.sur.pointData.class~=1), idsPlus);
ids = [idsPlus; idsMinus];
val = [ones(size(idsPlus)); -ones(size(idsMinus))];
o.lv.lap3 = solveLaplace(o.lv.L, ids, val, o.cfg.tol, o.cfg.maxit);

zeroCross3.cells = o.lv.vol.cells(any(diff(o.lv.lap3(o.lv.vol.cells)<0,1,2)~=0,2),:);
zeroCross3.cellTypes = repmat(uint8(10), size(zeroCross3.cells,1), 1);
pointIds = unique(zeroCross3.cells(:));
zeroCross3.points = o.lv.vol.points(pointIds,:);
zeroCross3.pointData.ids = pointIds;
zeroCross3.pointData.abLaplace = o.bv.abLaplace(o.lv.ids(pointIds));
zeroCross3.pointData.lap2 = o.lv.lap2(pointIds);
zeroCross3.pointData.lap3 = o.lv.lap3(pointIds);
zeroCross3.cells = changem(zeroCross3.cells, 1:size(zeroCross3.points,1), pointIds);

minusPi25 = vtkThreshold(vtkThreshold(vtkDataSetSurfaceFilter(zeroCross3), 'points', 'lap3', [-inf 0]), 'points', 'lap2', [-inf 0]);
plusPi25 = vtkThreshold(vtkThreshold(vtkDataSetSurfaceFilter(zeroCross3), 'points', 'lap3', [-inf 0]), 'points', 'lap2', [0 inf]);
idsMinusPi25 = unique([minusPi25.pointData.ids; vtkMapPointIds(o.lv.vol, juncEndoPost)]);
idsPlusPi25 = unique([plusPi25.pointData.ids; vtkMapPointIds(o.lv.vol, juncEndoAnt)]);

if o.cfg.exportLevel > 1
    vtkWrite(vtkDeleteDataArrays(minusPi25), sprintf('%slvMinusPi25.vtp', o.cfg.outPrefix));
    vtkWrite(vtkDeleteDataArrays(plusPi25), sprintf('%slvPlusPi25.vtp', o.cfg.outPrefix));
end

meanEpiAnt  = mean(o.lv.lap2(vtkMapPointIds(o.lv.vol, juncEpiAnt)));
meanEpiPost = mean(o.lv.lap2(vtkMapPointIds(o.lv.vol, juncEpiPost)));

idsPlus = unique([ ...
    intersect(o.lv.surToVol(o.lv.sur.pointData.class==3), find(o.lv.lap2>meanEpiPost & o.lv.lap2<meanEpiAnt)); ... % endo
    o.lv.surToVol(o.lv.sur.pointData.class==6); o.lv.surToVol(unique(o.lv.sur.cells(o.lv.sur.cellData.class==0,:))) ... % epi
    ]);
idsMinus = setdiff(o.lv.surToVol(o.lv.sur.pointData.class~=1), idsPlus);
ids = [idsPlus; idsMinus];
val = [ones(size(idsPlus)); -ones(size(idsMinus))];
o.lv.lap4 = solveLaplace(o.lv.L, ids, val, o.cfg.tol, o.cfg.maxit);

zeroCross4.cells = o.lv.vol.cells(any(diff(o.lv.lap4(o.lv.vol.cells)<0,1,2)~=0,2),:);
zeroCross4.cellTypes = repmat(uint8(10), size(zeroCross4.cells,1), 1);
pointIds = unique(zeroCross4.cells(:));
zeroCross4.points = o.lv.vol.points(pointIds,:);
zeroCross4.pointData.ids = pointIds;
zeroCross4.pointData.abLaplace = o.bv.abLaplace(o.lv.ids(pointIds));
zeroCross4.pointData.lap2 = o.lv.lap2(pointIds);
zeroCross4.pointData.lap4 = o.lv.lap4(pointIds);
zeroCross4.cells = changem(zeroCross4.cells, 1:size(zeroCross4.points,1), pointIds);

minusPi2 = vtkThreshold(vtkThreshold(vtkDataSetSurfaceFilter(zeroCross4), 'points', 'lap4', [0 inf]), 'points', 'lap2', [-inf 0]);
plusPi2 = vtkThreshold(vtkThreshold(vtkDataSetSurfaceFilter(zeroCross4), 'points', 'lap4', [0 inf]), 'points', 'lap2', [0 inf]);
idsPlusPi2 = plusPi2.pointData.ids;
idsMinusPi2 = minusPi2.pointData.ids;

if o.cfg.exportLevel > 1
    vtkWrite(vtkDeleteDataArrays(minusPi2), sprintf('%slvMinusPi2.vtp', o.cfg.outPrefix));
    vtkWrite(vtkDeleteDataArrays(plusPi2), sprintf('%slvPlusPi2.vtp', o.cfg.outPrefix));
end

%% final Laplace solution in the LV

idsZero     = setdiff(idsZero, [idsPlusPi25; idsPlusPi2; idsPlusPi; idsMinusPi; idsMinusPi2; idsMinusPi25]);
idsPlusPi25 = setdiff(idsPlusPi25, [idsPlusPi2; idsPlusPi; idsMinusPi; idsMinusPi2; idsMinusPi25]);
idsPlusPi2  = setdiff(idsPlusPi2, [idsPlusPi; idsMinusPi; idsMinusPi2; idsMinusPi25]);
idsPlusPi   = setdiff(idsPlusPi, [idsMinusPi; idsMinusPi2; idsMinusPi25]);
idsMinusPi  = setdiff(idsMinusPi, [idsMinusPi2; idsMinusPi25]);
idsMinusPi2 = setdiff(idsMinusPi2, idsMinusPi25);

ids = [idsZero; idsPlusPi25; idsPlusPi2; idsPlusPi; idsMinusPi; idsMinusPi2; idsMinusPi25];
val = [zeros(size(idsZero)); repmat(pi/2.5,size(idsPlusPi25)); repmat(pi/2,size(idsPlusPi2)); repmat(pi,size(idsPlusPi)); repmat(-pi,size(idsMinusPi)); repmat(-pi/2,size(idsMinusPi2)); repmat(-pi/2.5,size(idsMinusPi25))];
o.lv.rt = solveLaplace(o.lv.L, ids, val, o.cfg.tol, o.cfg.maxit);

if o.cfg.exportLevel > 1
    tmp = o.lv.vol;
    tmp.pointData.lap1 = single(o.lv.lap1);
    tmp.pointData.lap2 = single(o.lv.lap2);
    tmp.pointData.lap3 = single(o.lv.lap3);
    tmp.pointData.lap4 = single(o.lv.lap4);
    tmp.pointData.rt = single(o.lv.rt);
    vtkWrite(tmp, sprintf('%slv.vtu', o.cfg.outPrefix));
end

%% final Laplace solution in the RV

rvBound = vtkThreshold(o.rv.sur, 'cells', 'class', [0 0]);
rvBound.pointData.lap2 = o.lv.lap2(vtkMapPointIds(o.lv.vol, rvBound));
rvPlusPi2 = vtkThreshold(rvBound, 'points', 'lap2', [0 inf]);
rvMinusPi2 = vtkThreshold(rvBound, 'points', 'lap2', [-inf 0]);

idsRvPlusPi2 = vtkMapPointIds(o.rv.vol, rvPlusPi2);
idsRvMinusPi2 = vtkMapPointIds(o.rv.vol, rvMinusPi2);
idsToDelete = intersect(idsRvPlusPi2, idsRvMinusPi2);
idsRvPlusPi2 = setdiff(idsRvPlusPi2, idsToDelete);
idsRvMinusPi2 = setdiff(idsRvMinusPi2, idsToDelete);

if o.cfg.exportLevel > 1
    vtkWrite(vtkDeleteDataArrays(rvPlusPi2), sprintf('%srvPlusPi2.vtp', o.cfg.outPrefix));
    vtkWrite(vtkDeleteDataArrays(rvMinusPi2), sprintf('%srvMinusPi2.vtp', o.cfg.outPrefix));
end

ids = [idsRvPlusPi2; idsRvMinusPi2];
val = [repmat(pi/2,size(idsRvPlusPi2)); repmat(-pi/2,size(idsRvMinusPi2))];
o.rv.rt = solveLaplace(o.rv.L, ids, val, o.cfg.tol, o.cfg.maxit);

if o.cfg.exportLevel > 2
    tmp = o.rv.vol;
    tmp.pointData.rt = single(o.rv.rt);
    vtkWrite(tmp, sprintf('%srv.vtu', o.cfg.outPrefix));
end

%%

o.bv.rt = NaN(size(o.bv.vol.points,1),1);
o.bv.rt(o.lv.ids) = o.lv.rt;
o.bv.rt(o.rv.ids) = o.rv.rt;
bvNan = isnan(o.bv.rt);
lvNan = bvNan & o.bv.tv==0;
rvNan = bvNan & o.bv.tv==1;
lvNotNan = find(~bvNan & o.bv.tv==0);
rvNotNan = find(~bvNan & o.bv.tv==1);
o.bv.rt(lvNan) = o.bv.rt(lvNotNan(knnsearch(o.bv.vol.points(lvNotNan,:), o.bv.vol.points(lvNan,:))));
o.bv.rt(rvNan) = o.bv.rt(rvNotNan(knnsearch(o.bv.vol.points(rvNotNan,:), o.bv.vol.points(rvNan,:))));
o.bv.rtSin = sin(o.bv.rt);
o.bv.rtCos = cos(o.bv.rt);

o.result.pointData.rt = single(o.bv.rt);
o.result.pointData.rtSin = single(o.bv.rtSin);
o.result.pointData.rtCos = single(o.bv.rtCos);

if o.cfg.exportLevel > 0
    vtkWrite(o.result, sprintf('%sresult.vtu', o.cfg.outPrefix));
end

o.printStatus(sprintf('%.1f seconds\n', toc-t), true);
o.available.rotational = true;

end