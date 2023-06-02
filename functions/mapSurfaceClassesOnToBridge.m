function [struct1] = mapSurfaceClassesOnToBridge(struct1, struct2, struct3, struct4, struct5, cfg, mappingTol)
    % Map given boundary surfaces to surface of volume mesh and annotate classes
    %
    % [struct1] = mapSurfaceClassesOnToBridge(struct1, struct2, struct3, struct4, struct5, mappingTol)
    %
    % Inputs:
    %
    %   struct1, struct, tetrahedral vol mesh, with [numVolPoints x 3]
    %   struct2, struct, triangular sur mesh, with [Points x 3]
    %   struct3, struct, triangular sur mesh, with [Points x 3]
    %   struct4, struct, triangular sur mesh, with [Points x 3]
    %   struct5, struct, triangular sur mesh, with [Points x 3]
    %   cfg configuration prefix for saving the files produced, property of class cobiveco
    %   mappingTol, float,  mapping tolerance applied
    %
    % Output:
    %   
    %   struct1, struct, tetrahedral vol mesh, with [numVolPoints x 3]
    %
    % Written by Lisa Pankewitz

    % extract surface
    struct1.sur = vtkDataSetSurfaceFilter(struct1.vol);
    struct1.meanEdgLen = mean(vtkEdgeLengths(struct1.vol));   
    % apply surface to volume
    % mapping from surface to volume point ids
    struct1.surToVol = vtkMapPointIds(struct1.vol, struct1.sur); 
    boundaryStruct1SurAppended = vtkAppendPolyData({struct2,struct3, struct4, struct5});
    cnp = cumsum([size(struct2.points,1); size(struct3.points,1);size(struct4.points,1); size(struct5.points,1)]);
    ids = vtkMapPointIds(boundaryStruct1SurAppended, struct1.sur);
    struct1.sur.pointData.class = repmat(uint8(0),size(struct1.sur.points,1),1);
    struct1.sur.pointData.class(ids<=cnp(1)) = 1;               % struct2
    struct1.sur.pointData.class(ids>cnp(1) & ids<=cnp(2)) = 2;  % struct3
    struct1.sur.pointData.class(ids>cnp(2) & ids<=cnp(3)) = 3;  % struct4
    struct1.sur.pointData.class(ids>cnp(3)) = 4;  % struct5
    mappingDist = sqrt(sum((struct1.sur.points-boundaryStruct1SurAppended.points(ids,:)).^2,2));
    struct1.sur.pointData.class(mappingDist > mappingTol*struct1.meanEdgLen) = 0;
    
    struct1.vol.pointData.surClass = repmat(uint8(0),size(struct1.vol.points,1),1);
    struct1.vol.pointData.surClass(struct1.surToVol) = struct1.sur.pointData.class;
    vtkWrite(struct1.vol,sprintf('%sstruct1SurfaceClasses.vtu', cfg));


end