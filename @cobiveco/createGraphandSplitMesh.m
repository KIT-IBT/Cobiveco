function createGraphandSplitMesh(o)

% Creates a graph from the surface mesh of the object.
% Calcaulates the cutting planes by creating two shortest path
% between the inflow and outflow valve annuli  in each ventricle on the endo- and
% epicardial surface graphs and fitting planes through them.
%
% createGraphandSplitMesh(o)
%
% Input:
%   o: instance of the class cobiveco
%
% Outputs:
%   o.m0RvBridge, object of class cobiveco [struct] (For details see cobiveco class documentation)
%   o.m0LvBridge, object of class cobiveco [struct] (For details see cobiveco class documentation)
%   o.m0Ventricles, object of class cobiveco [struct] (For details see cobiveco class documentation)


if ~o.available.heartAxesAndApex
    o.computeHeartAxesAndApex;
end

o.printStatus('Creating Graph...');
t = toc;

%% CreateGraphAndSplitMesh
% extract edges
edgeStructM0 = vtkExtractEdges(o.m0.sur);

vtkWrite(edgeStructM0, sprintf('%sedgestructM0.vtp', o.cfg.outPrefix));

% remesh for defining septum
isovalue = 0.5;
numTries = 5;
for i = 1:numTries
    [m1remesh.vol,mmgStatus,mmgOutput2] = mmg(o.m1.vol, o.m1.ridgeLaplace, sprintf('-ls %1.5e -nr -hausd %1.5e -hmin %1.5e -hmax %1.5e', isovalue, o.cfg.mmgSizingParam(:)'*o.m0.meanEdgLen));

    if o.cfg.exportLevel > 1 || mmgStatus ~= 0
        if i == 1
            fid = fopen(sprintf('%smmgOutput_m1remesh.txt', o.cfg.outPrefix), 'w');
        else
            fid = fopen(sprintf('%smmgOutput_m1remesh.txt', o.cfg.outPrefix), 'a');
        end
        fprintf(fid, '%s', mmgOutput2);
        fclose(fid);
    end

    if mmgStatus == 0
        break;
    elseif i < numTries
        warning('Mmg remeshing with an isovalue of %.3f failed (system command status %i). Trying again with a slightly larger isovalue.', isovalue, mmgStatus);
        isovalue = isovalue + 0.1*prctile(max(o.m1.ridgeLaplace(o.m1.vol.cells),[],2)-min(o.m1.ridgeLaplace(o.m1.vol.cells),[],2),95);
        o.cfg.mmgSizingParam(1) = 0.8*o.cfg.mmgSizingParam(1);
    else
        error('Mmg remeshing failed ultimately after trying %i different isovalues (system command status %i).', numTries, mmgStatus);
    end
end

m1remesh.sur = vtkDataSetSurfaceFilter(vtkDeleteDataArrays(m1remesh.vol));
m1remesh.sur = vtkArrayMapperNearestNeighbor(o.m0.sur, m1remesh.sur);
m1remesh.surToVol = vtkMapPointIds(m1remesh.vol, m1remesh.sur);
vtkWrite(m1remesh.vol, sprintf('%sm1remeshVolCheckSeptum.vtu', o.cfg.outPrefix));


% Map celldata to point data
test = vtkCellDataToPointData(m1remesh.vol, (true));
mapClassValuesOntoM0 = vtkArrayMapperNearestNeighbor(test, o.m0.vol);
% find septal ids
idsSept = find(mapClassValuesOntoM0.pointData.class == 2);

% Create Graph for Surface Mesh
graphMesh0 = createGraphFromStruct(o.m0);

% Create Indexing array for SurToVol for Graph Nodes
graphMesh0EpicardialSurfaceOnlyNodeIndexArray = (1:1:size(o.m0.sur.points,1))';
volNodesOnM0correspondingToSurNodes = o.m0.surToVol(graphMesh0EpicardialSurfaceOnlyNodeIndexArray);

% Define indices of LV and RV Nodes including the Septum to be excluded in
% the (Sub)Graph Search
[idxRvInSurfaceNodes, idxLvInSurfaceNodes, oppositeIdxRvInSurfaceNodes, oppositeIdxLvInSurfaceNodes] = filterGlobalIndicesEndo(o.m0, volNodesOnM0correspondingToSurNodes,idsSept);

if o.cfg.exportLevel > 1
    checkExcludedNodesRv = zeros(size(o.m0.vol.points,1),1);
    checkExcludedNodesRv(o.m0.surToVol(idxRvInSurfaceNodes)) = 1;
    o.m0.vol.pointData.checkExcludedNodesRv = checkExcludedNodesRv;
    vtkWrite(o.m0.vol, sprintf('%scheckexcludednodes.vtu', o.cfg.outPrefix));
end


% find center of all valves
centerPv = vtkCenterOfArea(o.surPv);
centerTv = vtkCenterOfArea(o.surTv);
centerMv = vtkCenterOfArea(o.surMv);
centerAv = vtkCenterOfArea(o.surAv);

[tvPvVector, pvTvVector, mvAvVector, avMvVector] = createValveCenterConnectingVectors(centerMv, centerAv, centerTv,centerPv);

% calculate surface normals using pca
normalSurAv = pcaSurfaceNormal(o.surAv, o.longAx);
normalSurMv = pcaSurfaceNormal(o.surMv, o.longAx);
normalSurPv = pcaSurfaceNormal(o.surPv, o.longAx);
normalSurTv = pcaSurfaceNormal(o.surTv, o.longAx);

% MV: calculate the vector and minimizing the angle to mvAvVector
[referenceVectorMv, referenceVectorProjectedPlaneMv, referenceVectorArrayMv, idxSurPointMv] = calculateValveReferenceVector_update(o.surMv, normalSurAv, centerAv);

% AV: calculate the vector and minimizing the angle to avMvVector
[referenceVectorAv, referenceVectorProjectedPlaneAv, referenceVectorArrayAv, idxSurPointAv] = calculateValveReferenceVector_update(o.surAv, normalSurMv, centerMv);

% TV: calculate the vector and minimizing the angle to pvTvVector
[referenceVectorTv, referenceVectorProjectedPlaneTv, referenceVectorArrayTv, idxSurPointTv] = calculateValveReferenceVector_update(o.surTv,normalSurPv, centerPv);

% PV: calculate the vector and minimizing the angle to pvTvVector
[referenceVectorPv, referenceVectorProjectedPlanePv, referenceVectorArrayPv, idxSurPointPv] = calculateValveReferenceVector_update(o.surPv,normalSurTv,centerTv);

% project points on each valve annulus
[projectedPointsMv, referenceVectorArrayMv] = projectSurOnToPlane(o.surMv, normalSurMv);
[projectedPointsAv, referenceVectorArrayAv] = projectSurOnToPlane(o.surAv, normalSurAv);
[projectedPointsTv, referenceVectorArrayTv] = projectSurOnToPlane(o.surTv, normalSurTv);
[projectedPointsPv, referenceVectorArrayPv] = projectSurOnToPlane(o.surPv, normalSurPv);

% find the same point in the new projection for each valve based on the id
referenceVectorProjectedPlaneMv = referenceVectorArrayMv(idxSurPointMv,:);
referenceVectorProjectedPlaneAv = referenceVectorArrayAv(idxSurPointAv,:);
referenceVectorProjectedPlaneTv = referenceVectorArrayTv(idxSurPointTv,:);
referenceVectorProjectedPlanePv = referenceVectorArrayPv(idxSurPointPv,:);


% rotate vector around valve normal
[referenceVectorNewMvClockwise, R, t] = AxelRot(referenceVectorProjectedPlaneMv', 70, normalSurMv,[]);
[referenceVectorNewMvCounterclockwise, R, t] = AxelRot(referenceVectorProjectedPlaneMv', -70, normalSurMv,[]);
[referenceVectorNewAvClockwise, R, t] = AxelRot(referenceVectorProjectedPlaneAv', 70, normalSurAv,[]);
[referenceVectorNewAvCounterclockwise, R, t] = AxelRot(referenceVectorProjectedPlaneAv', -70, normalSurAv,[]);
[referenceVectorNewTvClockwise, R, t] = AxelRot(referenceVectorProjectedPlaneTv', 75, normalSurTv,[]);
[referenceVectorNewTvCounterclockwise, R, t] = AxelRot(referenceVectorProjectedPlaneTv', -20, normalSurTv,[]);
[referenceVectorNewPvClockwise, R, t] = AxelRot(referenceVectorProjectedPlanePv', 70, normalSurPv,[]);
[referenceVectorNewPvCounterclockwise, R, t] = AxelRot(referenceVectorProjectedPlanePv', -70, normalSurPv,[]);

% plot for axis
plotValveAxis(o.surMv,centerMv,normalSurMv,referenceVectorMv,referenceVectorNewMvClockwise',referenceVectorNewMvCounterclockwise', o.cfg.outPrefix, "Mv");
plotValveAxis(o.surAv,centerAv,normalSurAv,referenceVectorAv,referenceVectorNewAvClockwise',referenceVectorNewAvCounterclockwise', o.cfg.outPrefix, "Av");
plotValveAxis(o.surTv,centerTv,normalSurTv,referenceVectorTv,referenceVectorNewTvClockwise',referenceVectorNewTvCounterclockwise', o.cfg.outPrefix, "Tv");
plotValveAxis(o.surPv,centerPv,normalSurPv,referenceVectorPv,referenceVectorNewPvClockwise',referenceVectorNewPvCounterclockwise', o.cfg.outPrefix, "Pv");

% find vector best aligning (angle, distance) with the new (rotated) vectors
[pointMv80, pointMvMin80] = calculatePointsAligningWithRotatedVector(o.surMv, referenceVectorNewMvClockwise, referenceVectorNewMvCounterclockwise, referenceVectorArrayMv);
[pointAv80, pointAvMin80] = calculatePointsAligningWithRotatedVector(o.surAv, referenceVectorNewAvClockwise, referenceVectorNewAvCounterclockwise, referenceVectorArrayAv);
[pointTv80, pointTvMin80] = calculatePointsAligningWithRotatedVector(o.surTv, referenceVectorNewTvClockwise, referenceVectorNewTvCounterclockwise, referenceVectorArrayTv);
[pointPv80, pointPvMin80] = calculatePointsAligningWithRotatedVector(o.surPv, referenceVectorNewPvClockwise, referenceVectorNewPvCounterclockwise, referenceVectorArrayPv);

% Visualize starting points
visualizeStartingPointsPath(o.m0, pointAv80, pointAvMin80, pointMv80, pointMvMin80, pointTv80, pointTvMin80, pointPv80, pointPvMin80, o.cfg);

% Return Global indices of points in graph (could be exchanged for o.m0.sur)
[nodeIndexGlobalMv80, nodeIndexGlobalMvMin80] = findPointInGraph(edgeStructM0, pointMv80, pointMvMin80, oppositeIdxLvInSurfaceNodes);
[nodeIndexGlobalAv80, nodeIndexGlobalAvMin80] = findPointInGraph(edgeStructM0, pointAv80, pointAvMin80, oppositeIdxLvInSurfaceNodes);
[nodeIndexGlobalTv80, nodeIndexGlobalTvMin80] = findPointInGraph(edgeStructM0, pointTv80, pointTvMin80, oppositeIdxRvInSurfaceNodes);
[nodeIndexGlobalPv80, nodeIndexGlobalPvMin80] = findPointInGraph(edgeStructM0, pointPv80, pointPvMin80, oppositeIdxRvInSurfaceNodes);

% Create subgraph from graph containing only nodes of interest (rv, lv
% without septum/epi)
GsubRv = subgraph(graphMesh0, idxRvInSurfaceNodes);
GsubLv = subgraph(graphMesh0, idxLvInSurfaceNodes);

% map global node indices to local ones for each subgraph
nodeIndexLocalMvMin80 = mapGlobalIndexToGraphIndex(GsubLv, nodeIndexGlobalMvMin80, 'mvMin80');
nodeIndexLocalMv80 = mapGlobalIndexToGraphIndex(GsubLv, nodeIndexGlobalMv80, 'mv80');
nodeIndexLocalAvMin80 = mapGlobalIndexToGraphIndex(GsubLv, nodeIndexGlobalAvMin80, 'avMin80');
nodeIndexLocalAv80 = mapGlobalIndexToGraphIndex(GsubLv, nodeIndexGlobalAv80, 'av80');
nodeIndexLocalTvMin80 = mapGlobalIndexToGraphIndex(GsubRv, nodeIndexGlobalTvMin80, 'tvMin80');
nodeIndexLocalTv80 = mapGlobalIndexToGraphIndex(GsubRv, nodeIndexGlobalTv80, 'tv80');
nodeIndexLocalPvMin80 = mapGlobalIndexToGraphIndex(GsubRv, nodeIndexGlobalPvMin80, 'pvMin80');
nodeIndexLocalPv80 = mapGlobalIndexToGraphIndex(GsubRv, nodeIndexGlobalPv80, 'pv80');

% find the shortest path between two points in graph
path1MvAv = shortestpath(GsubLv, nodeIndexLocalMvMin80, nodeIndexLocalAv80);
path2MvAv = shortestpath(GsubLv, nodeIndexLocalMv80, nodeIndexLocalAvMin80);
path1TvPv = shortestpath(GsubRv, nodeIndexLocalTv80, nodeIndexLocalPvMin80);
path2TvPv = shortestpath(GsubRv, nodeIndexLocalTvMin80, nodeIndexLocalPv80);

path1TvPv = path1TvPv';
path2TvPv = path2TvPv';
path1MvAv = path1MvAv';
path2MvAv = path2MvAv';

% visualize paths calculated
arrayPath1TvPv = cellfun(@str2num,table2array(GsubRv.Nodes(path1TvPv,:)));
arrayPath2TvPv = cellfun(@str2num,table2array(GsubRv.Nodes(path2TvPv,:)));
arrayPath1MvAv = cellfun(@str2num,table2array(GsubLv.Nodes(path1MvAv,:)));
arrayPath2MvAv = cellfun(@str2num,table2array(GsubLv.Nodes(path2MvAv,:)));

o.m0.vol.pointData.shortestpath = zeros(size(o.m0.vol.points,1),1);
o.m0.vol.pointData.shortestpath(o.m0.surToVol(arrayPath1TvPv)) = 1;
o.m0.vol.pointData.shortestpath(o.m0.surToVol(arrayPath2TvPv)) = 2;
o.m0.vol.pointData.shortestpath(o.m0.surToVol(arrayPath1MvAv)) = 3;
o.m0.vol.pointData.shortestpath(o.m0.surToVol(arrayPath2MvAv)) = 4;

if o.cfg.exportLevel > 1
    vtkWrite(o.m0.vol, sprintf('%stestshortestpathLatest.vtu', o.cfg.outPrefix));
end

% Node IDs Epicardium
nodesEpi = o.m0.surToVol(o.m0.sur.pointData.class == 1 | o.m0.sur.pointData.class == 2);

% find closest coordinates in epicardium
[closestPointsInEpicardiumArrayPath1TvPv] = findMatchingPointinEpi(o.m0, arrayPath1TvPv, nodesEpi);
[closestPointsInEpicardiumArrayPath2TvPv] = findMatchingPointinEpi(o.m0, arrayPath2TvPv, nodesEpi);
[closestPointsInEpicardiumArrayPath1MvAv] = findMatchingPointinEpi(o.m0, arrayPath1MvAv, nodesEpi);
[closestPointsInEpicardiumArrayPath2MvAv] = findMatchingPointinEpi(o.m0, arrayPath2MvAv, nodesEpi);

if o.cfg.exportLevel > 1
    o.m0.vol.pointData.shortestpathepi = zeros(size(o.m0.vol.points,1),1);
    o.m0.vol.pointData.shortestpathepi(closestPointsInEpicardiumArrayPath1TvPv) = 1;
    o.m0.vol.pointData.shortestpathepi(closestPointsInEpicardiumArrayPath2TvPv) = 2;
    o.m0.vol.pointData.shortestpathepi(closestPointsInEpicardiumArrayPath1MvAv) = 3;
    o.m0.vol.pointData.shortestpathepi(closestPointsInEpicardiumArrayPath2MvAv) = 4;
    vtkWrite(o.m0.vol, sprintf('%stestshortestpathLatestEpipints.vtu', o.cfg.outPrefix));
end

% create field RegionIdVentricle
o.m0.vol.cellData.RegionIdVentricle = zeros(size(o.m0.vol.cells,1),1);
% fit plane to epi and respective endo values defined as paths
[o.m0, volIdxNonBridge1] = defineBridgesShortestPath(o.m0,arrayPath1TvPv,arrayPath2TvPv,closestPointsInEpicardiumArrayPath1TvPv, closestPointsInEpicardiumArrayPath2TvPv, o.cfg);
[o.m0, volIdxNonBridge2] = defineBridgesShortestPath(o.m0,arrayPath1MvAv,arrayPath2MvAv,closestPointsInEpicardiumArrayPath1MvAv, closestPointsInEpicardiumArrayPath2MvAv, o.cfg, 'rv');
% visualize
o.m0.vol.pointData.idxNonBridge = zeros(size(o.m0.vol.points,1),1);
o.m0.vol.pointData.idxNonBridge(volIdxNonBridge1) = 1;
o.m0.vol.pointData.idxNonBridge(volIdxNonBridge2) = 2;
vtkWrite(o.m0.vol, sprintf('%sbridgeArray.vtu', o.cfg.outPrefix));

% Remesh and split mesh into bridges and ventricles
[o.m0Ventricles ,o.m0RvBridge, o.m0LvBridge] = SplitMesh(o.m0, o.cfg);

% interpolate tv to current volume
Matrix = baryInterpMat(o.m0.vol.points, o.m0.vol.cells, o.m0Ventricles.vol.points);
o.m0Ventricles.tv = Matrix * o.m0.tv;
sur = vtkDataSetSurfaceFilter(o.m0Ventricles.vol);

% find surface classes and points
o.m0Ventricles.sur= vtkDataSetSurfaceFilter(vtkDeleteDataArrays(o.m0Ventricles.vol));
o.m0Ventricles.sur = vtkArrayMapperNearestNeighbor(o.m0.sur, o.m0Ventricles.sur);
o.m0Ventricles.surToVol = vtkMapPointIds(o.m0Ventricles.vol, o.m0Ventricles.sur);
o.m0Ventricles.meanEdgLen = mean(vtkEdgeLengths(o.m0Ventricles.vol));

% create class on volume
sur.surClass = repmat(uint8(0),size(o.m0Ventricles.sur.points,1),1);
sur.surClass = o.m0Ventricles.sur.pointData.class;

% figure out what are the tv values at the surface
% get volumen ids
sur.tv = o.m0Ventricles.tv(o.m0Ventricles.surToVol);
vtkWrite(sur, sprintf('%stestSlassificationsurface.vtp', o.cfg.outPrefix));

if o.cfg.exportLevel > 1
    testInfo = sur;
    testInfo.pointData.surClass = o.m0Ventricles.sur.pointData.class;
    testInfo.pointData.tv =  sur.tv;
end

vtkWrite(testInfo, sprintf('%sdebug.vtp', o.cfg.outPrefix));

% extract only surface nodes
o.m0RvBridge.sur= vtkDataSetSurfaceFilter(vtkDeleteDataArrays(o.m0RvBridge.vol));
o.m0LvBridge.sur= vtkDataSetSurfaceFilter(vtkDeleteDataArrays(o.m0LvBridge.vol));
o.m0Ventricles.sur= vtkDataSetSurfaceFilter(vtkDeleteDataArrays(o.m0Ventricles.vol));

if o.cfg.exportLevel > 1
    vtkWrite(o.m0RvBridge.sur, sprintf('%srvbridgesurcheck.vtp', o.cfg.outPrefix));
    vtkWrite(o.m0LvBridge.sur,sprintf('%slvbridgesurcheck.vtp', o.cfg.outPrefix));
    vtkWrite(o.m0Ventricles.sur,sprintf('%sventricularsurcheck.vtp', o.cfg.outPrefix));
end

o.printStatus(sprintf('%.1f seconds\n', toc-t), true);
o.available.creategraphandsplitmesh = true;

end
