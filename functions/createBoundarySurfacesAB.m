function [struct1, struct1FwBoundary, struct1SeptumBoundary] = createBoundarySurfacesAB(struct1, struct2, referenceVectorFW, referenceVectorSeptum, longAx, cfg)
    % Create 2 boundary surface for the fsw coordinate that you can apply readily to the bridge.
    % checks which boundary belongs to which site of the bridge (free wall or septum), by comparing it to the reference vector.
    %
    % [struct1, struct1FwBoundary, struct1SeptumBoundary] = createBoundarySurfacesAB(struct1, struct2, referenceVectorFW, referenceVectorSeptum, cfg)
    %
    % Inputs:
    %
    %   struct1, struct, tetrahedral vol mesh, with [numVolPoints x 3]
    %   struct2, struct, tetrahedral vol mesh, with [numVolPoints x 3]
    %   referenceVectorFW [1x3]
    %   referenceVectorSeptum [1x3]
    %
    % Output:
    %
    %   struct1, struct, triangular sur mesh, with [Points x 3]
    %   struct2, struct, triangular sur mesh, with [Points x 3]
    %   cfg, object/structure in cobiveco class which contains info about configuration storage location
    %
    %
    % Written by Lisa Pankewitz
    % 
    condaPath = '/Users/lieschen/anaconda3/bin';    % replace this with the path to your conda installation, see README.md
    conda.addBaseCondaPath(condaPath)
    % extract surfaces
    struct1.sur= vtkDataSetSurfaceFilter(vtkDeleteDataArrays(struct1.vol));
    struct2.sur= vtkDataSetSurfaceFilter(vtkDeleteDataArrays(struct2.vol));

    %  check that the surface extracted is actually right
    vtkWrite(struct1.sur,sprintf('%sstruct1.vtp', cfg));
    vtkWrite(struct2.sur,sprintf('%sstruct2.vtp', cfg));

    % find the  initial boundary surfaces on struct 1
    struct1InitialBoundaryIndices = knnsearch(struct1.sur.points, struct2.sur.points, 'NSMethod','kdtree');
    struct1.sur.pointData.keep = repmat(uint8(0),size(struct1.sur.points,1),1);
    indicesStruct1ToKeep = unique(struct1InitialBoundaryIndices);
    struct1.sur.pointData.keep(indicesStruct1ToKeep) = 1;
    struct1InitialBoundary  = vtkDeleteDataArrays(vtkThreshold(struct1.sur, 'points', 'keep', [1 1]));
    vtkWrite(struct1InitialBoundary,sprintf('%sstruct1InitialBoundary.vtp', cfg));

    % find the  initial boundary surfaces on struct 2
    struct2InitialBoundaryIndices = knnsearch(struct2.sur.points,struct1.sur.points,  'NSMethod','kdtree');
    struct2.sur.pointData.keep = repmat(uint8(0),size(struct2.sur.points,1),1);
    IndicesStruct2ToKeep = unique(struct2InitialBoundaryIndices);
    struct2.sur.pointData.keep(IndicesStruct2ToKeep) = 1;
    struct2InitialBoundary  = vtkDeleteDataArrays(vtkThreshold(struct2.sur, 'points', 'keep', [1 1]));
    vtkWrite(struct2InitialBoundary,sprintf('%sstruct2InitialBoundary.vtp', cfg));


    %reiterate using initial boundaries
    struct1ReiteratedBoundaryIndices = knnsearch(struct1InitialBoundary.points, struct2InitialBoundary.points,'NSMethod','kdtree');
    struct1InitialBoundary.pointData.keep = repmat(uint8(0),size(struct1InitialBoundary.points,1),1);
    indicesStruct1ToKeep = unique(struct1ReiteratedBoundaryIndices);
    struct1InitialBoundary.pointData.keep(indicesStruct1ToKeep) = 1;
    struct1Boundary = vtkDeleteDataArrays(vtkThreshold(struct1InitialBoundary, 'points', 'keep', [1 1]));
    vtkWrite(struct1Boundary,sprintf('%sstruct1Boundary.vtp', cfg));

    % use connectivity filter to find the 2 biggest boundaries
    struct1BoundaryRegions = vtkConnectivityFilter(struct1Boundary);

    % find out which are the two biggest ones by looking at the group count
    [GroupCount, GroupRegion] = groupcounts(struct1BoundaryRegions.pointData.RegionId);

    if ~(isempty(GroupCount)) && size(GroupCount,1) >= 2
        disp('Found enough boundaries. Continue.');
    else
        error('Could not locate at least two boundaries. Only found %i. Please check manually.', size(GroupCount,1));
    end


    % order the group count
    [GroupCountValues, originalIndices] = sort(GroupCount,'descend');

    % Take the two biggest counts
    boundary1Index = int64(GroupRegion(originalIndices(1)));
    boundary1Index = cast(boundary1Index,'double');
    boundary2Index = GroupRegion(originalIndices(2));
    boundary2Index = cast(boundary2Index,'double');
    disp(boundary1Index)
    boundaryRegions1 = struct1BoundaryRegions;
    boundaryRegions2 = struct1BoundaryRegions;
    struct1Boundary1 = vtkDeleteDataArrays(vtkThreshold(boundaryRegions1, 'points', 'RegionId', [boundary1Index boundary1Index]));
    struct1Boundary2 = vtkDeleteDataArrays(vtkThreshold(boundaryRegions2, 'points', 'RegionId', [boundary2Index boundary2Index]));
    vtkWrite(struct1Boundary1,sprintf('%sstruct1Boundary1.vtp', cfg));
    vtkWrite(struct1Boundary2,sprintf('%sstruct1Boundary2.vtp', cfg));

    % now we check if the boundary is at the FW or septum
    meanBoundary1 = mean(struct1Boundary1.points);
    meanBoundary2 = mean(struct1Boundary2.points);
    % vector boundary one to two
    vectorB1B2 = meanBoundary2 - meanBoundary1;
    % vector boundary two to one
    vectorB2B1 = meanBoundary1 - meanBoundary2;

    if dot(vectorB1B2,referenceVectorSeptum) > 0 && dot(vectorB2B1,referenceVectorSeptum) < 0
        struct1BoundaryFw = struct1Boundary1;
        struct1BoundarySeptum = struct1Boundary2;

    elseif dot(vectorB1B2,referenceVectorSeptum) < 0 && dot(vectorB2B1,referenceVectorSeptum) > 0
        struct1BoundaryFw = struct1Boundary2;
        struct1BoundarySeptum = struct1Boundary1;

    else
        error('Cannot Match boundaries to FW and septum. Please check manually.')
    end

    vtkWrite(struct1BoundaryFw,sprintf('%sstruct1BoundaryFw.vtp', cfg));
    vtkWrite(struct1BoundarySeptum,sprintf('%sstruct1BoundarySeptum.vtp', cfg));

    % find these exact ones in the original bridge so we can apply them, 1 for fw and 2 for septum
    struct1BoundaryFwIndices = knnsearch(struct1.sur.points, struct1BoundaryFw.points, 'NSMethod','kdtree');
    struct1BoundarySeptumIndices = knnsearch(struct1.sur.points, struct1BoundarySeptum.points, 'NSMethod','kdtree');
    struct1.sur.pointData.boundaries = repmat(uint8(0),size(struct1.sur.points,1),1);
    indicesStruct1FwToKeep = unique(struct1BoundaryFwIndices);
    indicesStruct1SeptumToKeep = unique(struct1BoundarySeptumIndices);
    struct1.sur.pointData.boundaries(indicesStruct1FwToKeep) = 1;
    struct1.sur.pointData.boundaries(indicesStruct1SeptumToKeep) = 2;
    struct1FwBoundary  = vtkDeleteDataArrays(vtkThreshold(struct1.sur, 'points', 'boundaries', [1 1]));
    struct1SeptumBoundary  = vtkDeleteDataArrays(vtkThreshold(struct1.sur, 'points', 'boundaries', [2 2]));
    vtkWrite(struct1FwBoundary,sprintf('%sstruct1FwBoundary.vtp', cfg));
    vtkWrite(struct1SeptumBoundary,sprintf('%sstruct1SeptumBoundary.vtp', cfg));

    % if there are points over a certain threshold close to the midpoint
    % between boundary and middle of bridges, find index of the elements and remove them
    % using bfs

    % vector FW direction
    vectorTowardsFW = -vectorB1B2;
    % vector septal direction
    vectorTowardsS = vectorB1B2;
    % create matrices based on reference
    referenceMatrixFWDir = repmat(vectorTowardsFW,size(struct1SeptumBoundary.points,1),1);
    referenceMatrixSDir = repmat(vectorTowardsS,size(struct1FwBoundary.points,1),1);

    % find mean point points
    meanPointS = mean(struct1SeptumBoundary.points,1);
    meanPointFW = mean(struct1FwBoundary.points,1);
    [idx, member] =  ismembertol(struct1SeptumBoundary.points, meanPointS, 3,'ByRows',true,'DataScale',[1,1,1]);
    visualizeMean = find(idx==1);
    % visualize
    struct1SeptumBoundary.pointData.meanSeptalBoundary = zeros(size(struct1SeptumBoundary.points,1));
    struct1SeptumBoundary.pointData.meanSeptalBoundary(visualizeMean) = 1;
    vtkWrite(struct1SeptumBoundary,sprintf('%sstruct1SeptumBoundaryMean.vtp', cfg));

    [idx, member] =  ismembertol(struct1FwBoundary.points, meanPointFW, 3,'ByRows',true,'DataScale',[1,1,1]);
    visualizeMean = find(idx==1);
    % visualize
    struct1FwBoundary.pointData.meanFwBoundary = zeros(size(struct1FwBoundary.points,1));
    struct1FwBoundary.pointData.meanFwBoundary(visualizeMean) = 1;
    vtkWrite(struct1FwBoundary,sprintf('%sstruct1FwBoundaryMean.vtp', cfg));

    % Look at the distance
    % find point between the two bridges
    pointInBetweenBridges = (meanBoundary2 + meanBoundary1)/2;

    % pre filter
    % start with septal boundary
    referencePointMatrixS = repmat(meanPointS,size(struct1SeptumBoundary.points,1),1);
    % create vectors
    comparisonVectorsS = struct1SeptumBoundary.points - referencePointMatrixS;
    % determine how well vector aligns with reference
    dotProductResultS = dot(comparisonVectorsS,referenceMatrixFWDir,2);
    alignedPointsS = find(dotProductResultS <0);

    % visualize check
    struct1SeptumBoundary.pointData.guess = zeros(size(struct1SeptumBoundary.points,1));
    struct1SeptumBoundary.pointData.guess(alignedPointsS) = 1;
    vtkWrite(struct1SeptumBoundary,sprintf('%sstruct1SeptumBoundary.vtp', cfg));
    vtkWrite(struct1SeptumBoundary,sprintf('%sstruct1SeptumBoundary.vtk', cfg));

    % make array of mean
    referencePointMatrixFW = repmat(meanPointFW,size(struct1FwBoundary.points,1),1);
    % create vectors
    comparisonVectorsFW = struct1FwBoundary.points - referencePointMatrixFW;
    % determine how well vector aligns with reference
    dotProductResultFW = dot(comparisonVectorsFW,referenceMatrixSDir,2);
    alignedPointsFW = find(dotProductResultFW <0);

    %visualize check
    struct1FwBoundary.pointData.guess = zeros(size(struct1FwBoundary.points,1));
    struct1FwBoundary.pointData.guess(alignedPointsFW) = 1;
    vtkWrite(struct1FwBoundary,sprintf('%sstruct1FwBoundary.vtp', cfg));

    subsetPointsFw = struct1FwBoundary.points(alignedPointsFW,:);
    subsetPointsS = struct1SeptumBoundary.points(alignedPointsS,:);
    % compare distance to this point in the bridge boundaries
    distancesFwFromMiddlePoint = pdist2(subsetPointsFw,pointInBetweenBridges,'euclidean');
    distancesSFromMiddlePoint = pdist2(subsetPointsS,pointInBetweenBridges,'euclidean');

    distanceMeanFw = pdist2(meanPointFW,pointInBetweenBridges);
    distanceMeanS = pdist2(meanPointS,pointInBetweenBridges);

    % normalize
    normalizedDistanceFw = distancesFwFromMiddlePoint/distanceMeanFw;
    normalizedDistanceS = distancesSFromMiddlePoint/distanceMeanS;

    % check if we have extra element in the middle of the bridge that need
    % to be removed to improve the accuracy of the boundary
    % start with fw boundary

    if ~isempty(normalizedDistanceFw(normalizedDistanceFw<0.5))
        disp('Need to remove some extra triangles.');
        % minimal distance to this point
        idInSubsetPointsFw = find(distancesFwFromMiddlePoint == min(distancesFwFromMiddlePoint));
        startingPointFw = alignedPointsFW(idInSubsetPointsFw);

        % Visualize
        struct1FwBoundary.pointData.initial = zeros(size(struct1SeptumBoundary.points,1));
        struct1FwBoundary.pointData.initial(startingPointFw) = 1;
        vtkWrite(struct1FwBoundary,sprintf('%sstruct1FwBoundaryStartingpoint.vtk', cfg));
        % export as vtk
        commandMeshio = sprintf('meshio convert %sstruct1FwBoundaryUnstr.vtk %sstruct1FwBoundaryUnstr.ply --ascii', cfg, cfg);
        % call BFS program in python to define extra triangles
        commandPython = sprintf('python extractExtraTriangles.py %d %sstruct1FwBoundaryUnstr.ply %sstructRemovedFw ',startingPointFw,cfg,cfg);
        [status,cmdout] = system(commandMeshio)
        [status,cmdout] = system(commandPython)

        partsToBeRemoved = vtkRead('out.id1.ply');

        toBeRemoved = knnsearch(struct1FwBoundary.points, partsToBeRemoved.points, 'NSMethod','kdtree');
        struct1FwBoundary.pointData.remove = repmat(uint8(0),size(struct1FwBoundary.points,1),1);
        indicesStruct1ToKeep = unique(toBeRemoved);
        struct1FwBoundary.pointData.remove(indicesStruct1ToKeep) = 1;
        struct1InitialBoundaryFw = vtkDeleteDataArrays(vtkThreshold(struct1FwBoundary, 'points', 'remove', [0 0]));
        vtkWrite(struct1InitialBoundaryFw,sprintf('%sstruct1FinalBoundaryFw.vtp', cfg));
        struct1FwBoundary = struct1InitialBoundaryFw;

    else
        disp('Thanks for checking. No need to do anything here.');
    end

    if ~isempty(normalizedDistanceS(normalizedDistanceS<0.5))
        disp('Need to remove some extra triangles.');
        % minimal distance to this point
        idInSubsetPointsSeptum = find(distancesSFromMiddlePoint == min(distancesSFromMiddlePoint));
        startingPointSeptum = alignedPointsS(idInSubsetPointsSeptum);

        % Visualize
        struct1SeptumBoundary.pointData.initial = zeros(size(struct1SeptumBoundary.points,1));
        struct1SeptumBoundary.pointData.initial(startingPointSeptum) = 1;
        % polydata to unstr grid
        vtkWrite(struct1SeptumBoundary,sprintf('%sstruct1SeptumBoundaryStartingpoint.vtk', cfg));
        % deactivate conda and use pv python to convert polydata to unstructured grid
        conda.setenv('base')
        conda.setenv('ldrb')
        conda.deactivate()


        commandPvpython = sprintf('pvpython ./../functions/polydata2unstrucgrid.py struct1SeptumBoundary.vtk %sstruct1SeptumBoundary.vtk %sstruct1SeptumBoundaryUnstr.vtk',cfg, cfg)
        [status,cmdout] = system(commandPvpython)
        if status ~= 0
            error('Paraview did not successfully transform polydata to unstructured grid. Need manual handling. Exiting.');
        end
        % set conda env that has meshio already installed
        conda.setenv('base')
        % export as vtk
        commandMeshio = sprintf('meshio convert %sstruct1SeptumBoundaryUnstr.vtk %sstruct1SeptumBoundaryUnstr.ply --ascii', cfg, cfg);
        % call BFS program in python to define extra triangles
        commandPython = sprintf('python extractExtraTriangles.py %d %sstruct1SeptumBoundaryUnstr.ply %sstructRemovedSeptum ',startingPointSeptum,cfg,cfg);
        [status,cmdout] = system(commandMeshio)
        if status ~= 0
            error('Meshio did not succeccfully transform .vtk to .ply. Need manual handling. Exiting.');
        end
        [status,cmdout] = system(commandPython)
        if status ~= 0
            error('Error in determining surface to be removed in extractExtraTriangles.py. Exiting.');
        end

        partsToBeRemoved = vtkRead(sprintf('%sstructRemovedSeptum.ply',cfg));
        toBeRemoved = knnsearch(struct1SeptumBoundary.points, partsToBeRemoved.points, 'NSMethod','kdtree');
        struct1SeptumBoundary.pointData.remove = repmat(uint8(0),size(struct1SeptumBoundary.points,1),1);
        indicesStruct1ToKeep = unique(toBeRemoved);
        struct1SeptumBoundary.pointData.remove(indicesStruct1ToKeep) = 1;
        struct1InitialBoundarySeptum = vtkDeleteDataArrays(vtkThreshold(struct1SeptumBoundary, 'points', 'remove', [0 0]));
        vtkWrite(struct1InitialBoundarySeptum,sprintf('%sstruct1FinalBoundaryS.vtp', cfg));

        struct1SeptumBoundary = struct1InitialBoundarySeptum;
    else
        disp('Thanks for checking. No need to do anything here.');
    end



end
