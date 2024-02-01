function [mesh, pointIdsBridges] = defineBridgesShortestPath(mesh, path1, path2, epiPath1, epiPath2, cfg, tv)

    % Defines the shortest path for bridges.
    %
    % [mesh, pointIdsBridges] = defineBridgesShortestPath(mesh, path1, path2, epiPath1, epiPath2, cfg, tv)
    % 
    % Inputs:
    %
    %   struct [numSurPoints x 3]
    %   path1, array double [nx1] (surface nodes)
    %   path2, array double [nx1] (surface nodes)
    %   epiPath1, array double [nx1]
    %   epiPath2, array double [nx1]
    %   cfg, object/structure in cobiveco class 
    %       which contains info about configuration storage location
    %   tv, boolean value where true represents tv and false represents lv. 
    %       If not given as input argument, tv is set to false
    %
    % Output:
    %
    %   mesh 
    %   pointIdsBridges, double [nx1]
    %
    % Written by Lisa Pankewitz

    if nargin < 7
        tv =  false;
    end
    
    % Fit plane through points
    % find points in plane
    idxPointsPlane1 =  mesh.surToVol(path1);
    % add epi to assumptions for the plane
    idxPointsPlane1 =  [idxPointsPlane1;epiPath1];
    pointsPlane1 = mesh.vol.points(idxPointsPlane1,:);
    plane1 = fit([double(mesh.vol.points(idxPointsPlane1,1)), double(mesh.vol.points(idxPointsPlane1,2))],double(mesh.vol.points(idxPointsPlane1,3)),'cubicinterp');

    idxPointsPlane2 =  mesh.surToVol(path2);
    idxPointsPlane2 =  [idxPointsPlane2;epiPath2];
    pointsPlane2 = mesh.vol.points(idxPointsPlane2,:);
    plane2 = fit([double(mesh.vol.points(idxPointsPlane2,1)), double(mesh.vol.points(idxPointsPlane2,2))],double(mesh.vol.points(idxPointsPlane2,3)),'cubicinterp');

    crossproduct = cross(pointsPlane1(2,:) - pointsPlane1(1,:),pointsPlane1(3,:) - pointsPlane1(1,:));
    normalPlane1 = crossproduct/norm(crossproduct);
    pc1 = pca(pointsPlane1);
    normalPlane1 = pc1(1,:);

    crossproduct = cross(pointsPlane2(2,:) - pointsPlane2(1,:),pointsPlane2(3,:) - pointsPlane2(1,:));
    normalPlane2 = crossproduct/norm(crossproduct);
    pc2 = pca(pointsPlane2);
    normalPlane2 = pc2(1,:);


    % get direction fro a point in plane1 to plane2 to get the major direction
    vetorReferenceNormalPlane1 = mean(pointsPlane1) - mean(pointsPlane2);

    % check the dot product between the reference and the plane normal 1, if negative, flip
    controlPlane1 = dot(vetorReferenceNormalPlane1,normalPlane1);
    controlPlane2 = dot(vetorReferenceNormalPlane1,normalPlane2);


    if controlPlane1 < 0
        normalPlane1 = -normalPlane1;
    end

    if controlPlane2 >= 0
        normalPlane2 = -normalPlane2;
    end

    % now find all the points  that are between the 2 planes

    if cfg.exportLevel > 2

        surAxes.points = [pointsPlane1(1,:); pointsPlane1(1,:)+5*pc1(1,:); pointsPlane1(1,:)+5*pc1(3,:); pointsPlane1(2,:)+3*pc1(2,:)];
        surAxes.cells = int32([1 2; 1 3; 1 4]);
        surAxes.cellTypes = uint8([3; 3; 3]);
        surAxes.cellData.axis = uint8([1; 2; 3]);
        vtkWrite(surAxes, sprintf('%stestbridgeaxisdefintion.vtp', cfg.outPrefix));
    
    end

    x = mesh.vol.points(:,1);
    y = mesh.vol.points(:,2);

    zInplane1 = plane1(double(x), double(y));
    % concat points
    pointsInPlane1 = horzcat(x,y,zInplane1);
    zInplane2 = plane2(double(x), double(y));
    % concat points
    pointsInPlane2 = horzcat(x,y,zInplane2);

    % check points in plane if nan look for closest point and exchange index
    % for that one and add to direction calculation to figure out if the point
    % is above the plane or not
    % prefilter all rows that have at least one nan element
    pointsInPlane1Filtered = pointsInPlane1(all(~isnan(pointsInPlane1),2),:);
    pointsInPlane2Filtered = pointsInPlane2(all(~isnan(pointsInPlane2),2),:);


    % Find max distance between points in plane and points in plane 2
    [index, dist] = knnsearch(pointsInPlane1Filtered,pointsInPlane2Filtered);
    maxDistance = max(dist);

    [index1, dist1] = knnsearch(mesh.vol.points,pointsInPlane1Filtered);
    [index2, dist2] = knnsearch(mesh.vol.points,pointsInPlane2Filtered);

    % use tenths percentile as threshold for distance search
    threshold = max(dist2)+prctile(dist2,75);

    %[idx1,D1] = rangesearch(mesh.vol.points,pointsInPlane1Filtered,1.0,'NSMethod','exhaustive');
    [idx1,D1] = rangesearch(mesh.vol.points,pointsInPlane1Filtered,threshold,'NSMethod','exhaustive');
    tf = ~cellfun(@isempty, idx1);

	idx1Concat = [idx1{tf}] ;
    idx1Concat = idx1Concat';

    %[idx2,D2] = rangesearch(mesh.vol.points,pointsInPlane2Filtered,1.0,'NSMethod','exhaustive');
    [idx2,D2] = rangesearch(mesh.vol.points,pointsInPlane2Filtered,threshold,'NSMethod','exhaustive');
    tf2 = ~cellfun(@isempty, idx2);

	idx2Concat = [idx2{tf2}] ;
    idx2Concat = idx2Concat';


    testIsMemberTolPlane2 = ismembertol(mesh.vol.points,pointsInPlane2Filtered,1.5,'ByRows',true,'DataScale',1.0);
    [indicesPlane2,values] = find(testIsMemberTolPlane2==1);

    closestPointsInPlane1 = mesh.vol.points(idx1Concat,:);
    closestPointsInPlane2 = mesh.vol.points(idx2Concat,:);

    % exclude all rv  or lv points
    closestPoints = zeros(size(mesh.vol.points,1),1);
    mesh.vol.pointData.closestpoints1 = closestPoints;
    mesh.vol.pointData.closestpoints1(idx2Concat) = 2;
    mesh.vol.pointData.closestpoints1(idx1Concat) = 1;
    vtkWrite(mesh.vol, sprintf('%sclosepoints.vtu', cfg.outPrefix));

    % exclude rv or lv points respectively
    if tv == false
        idsTv = find(mesh.tv == 0);
        updatedClosesPoint1 = setdiff(idx1Concat,idsTv);
        updatedClosesPoint2 = setdiff(idx2Concat,idsTv);
    else
        idsTv = find(mesh.tv == 1);
        updatedClosesPoint1 = setdiff(idx1Concat,idsTv);
        updatedClosesPoint2 = setdiff(idx2Concat,idsTv);
    end

    % exclude points from mesh
    closestPoints = zeros(size(mesh.vol.points,1),1);
    mesh.vol.pointData.closestpoints1 = closestPoints;
    mesh.vol.pointData.closestpoints1(updatedClosesPoint1) = 1;
    mesh.vol.pointData.closestpoints1(updatedClosesPoint2) = 1;
    vtkWrite(mesh.vol, sprintf('%sclosepointsFiltered.vtu', cfg.outPrefix));

    % created filter for thresholding
    bridge = mesh.vol;
    bridge.pointData.connector = int16(mesh.vol.pointData.closestpoints1);
    bridgeTest = vtkDeleteDataArrays(vtkThreshold(bridge, 'points', 'connector', [0 0]));
    vtkWrite(bridgeTest, sprintf('%sclosepointsFilteredBridgetest.vtu', cfg.outPrefix));

    % find the surfaces connected
    bridgeTestRegions = vtkConnectivityFilter(bridgeTest);
    [GC, GR] = groupcounts(bridgeTestRegions.pointData.RegionId);


    % check which region is bigger using the group count, which has the
    % same indexing as the group regions
    if GC(1) > GC(2)
        ventricle = vtkDeleteDataArrays(vtkThreshold(bridgeTestRegions, 'points', 'RegionId', [0 0]));
        bridgeTestRegion = vtkDeleteDataArrays(vtkThreshold(bridgeTestRegions, 'points', 'RegionId', [1 1]));
    else
        ventricle = vtkDeleteDataArrays(vtkThreshold(bridgeTestRegions, 'points', 'RegionId', [1 1]));
        bridgeTestRegion = vtkDeleteDataArrays(vtkThreshold(bridgeTestRegions, 'points', 'RegionId', [0 0]));
    end

    vtkWrite(bridgeTestRegion, sprintf('%sbridgeFiltered.vtu', cfg.outPrefix));
    vtkWrite(ventricle, sprintf('%sventricle.vtu', cfg.outPrefix));

    idsVentricle = find(bridgeTestRegions.pointData.RegionId == 1);
    bridgeTestRegions.pointData.RegionId(idsVentricle) = 10 ;

    mapVentricleRegion = vtkArrayMapperNearestNeighbor(bridgeTestRegions, mesh.vol);
    if GC(1) > GC(2)

        idsNotVentricle2 = find(mapVentricleRegion.pointData.RegionId ~=0);
        idsNotVentricle1 = find(bridge.pointData.connector ~=0);

    else
        idsNotVentricle2 = find(mapVentricleRegion.pointData.RegionId ==0);
        idsNotVentricle1 = find(bridge.pointData.connector ~=0);

    end

    idsNotVentricle = [idsNotVentricle1;idsNotVentricle2];

    if eval('isfield(mesh.vol.pointData,''RegionIdVentricle'')','0')  && ~isempty(mesh.vol.pointData.RegionIdVentricle)
        disp('All is well. RegionIdVentricle does exist.');
    else
        mesh.vol.pointData.RegionIdVentricle = zeros(size(mesh.vol.points,1),1);
    end

    if tv == false
         mesh.vol.pointData.RegionIdVentricle(idsNotVentricle) = 1;
    else
        mesh.vol.pointData.RegionIdVentricle(idsNotVentricle) = 2;
    end

    elementsSharingNodeWithVentricle = find(any(ismember(mesh.vol.cells,idsNotVentricle),2));

    if eval('isfield(mesh.vol.cellData,''RegionIdVentricle'')','0')  && ~isempty(mesh.vol.cellData.RegionIdVentricle)
        disp('All is well. RegionIdVentricle does exist.');
    else
        mesh.vol.cellData.RegionIdVentricle = zeros(size(mesh.vol.cells,1),1);
    end

    % rv will be set to 1
    % lv will be set to 2
    if tv == false
        mesh.vol.cellData.RegionIdVentricle(elementsSharingNodeWithVentricle) = 1;
        bridgeTestComple = vtkDeleteDataArrays(vtkThreshold(mesh.vol, 'cells', 'RegionIdVentricle', [1 1]));
        bridgeTestComple2 = vtkDeleteDataArrays(vtkThreshold(mesh.vol, 'cells', 'RegionIdVentricle', [0 0]));
    else
        mesh.vol.cellData.RegionIdVentricle(elementsSharingNodeWithVentricle) = 2;
        bridgeTestComple = vtkDeleteDataArrays(vtkThreshold(mesh.vol, 'cells', 'RegionIdVentricle', [2 2]));
        bridgeTestComple2 = vtkDeleteDataArrays(vtkThreshold(mesh.vol, 'cells', 'RegionIdVentricle', [0 0]));
    end

    % cell2pointdata
    outStruct = vtkCellDataToPointData(mesh.vol, (true));
    pointIdsBridges = find(outStruct.pointData.RegionIdVentricle == 1);

end
