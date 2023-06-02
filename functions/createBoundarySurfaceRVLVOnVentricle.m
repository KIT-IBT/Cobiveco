function [ventricles, RVBoundary, LVBoundary] = createBoundarySurfaceRVLVOnVentricle(ventricles, RVBoundaryFW, RVBoundarySeptum, LVBoundaryFW, LVBoundarySeptum, cfg)
    % Create boundary for RV and LV on ventricular volume
    %
    % [ventricles, RVBoundary, LVBoundary] = createBoundarySurfaceRVLVOnVentricle(ventricles, RVBoundaryFW, RVBoundarySeptum, LVBoundaryFW, LVBoundarySeptum, cfg)
    %
    % Inputs:
    %
    %   ventricles,         struct, tetrahedral vol mesh, with [numVolPoints x 4]
    %   RVBoundaryFW,       struct, triangular sur mesh, with [Points x 3]
    %   RVBoundarySeptum,   struct, triangular sur mesh, with [Points x 3]
    %   LVBoundaryFW,       struct, triangular sur mesh, with [Points x 3]
    %   LVBoundarySeptum,   struct, triangular sur mesh, with [Points x 3]
    %   cfg configuration prefix for saving the files produced
    %
    % Output:
    %
    %   ventricles,         struct, tetrahedral vol mesh, with [numVolPoints x 4]
    %   RVBoundary, struct, triangular sur mesh, with [Points x 3]
    %   LVBoundary, struct, triangular sur mesh, with [Points x 3]
    %
    %
    % Written by Lisa Pankewitz

    % extract surfaces
    ventricles.sur= vtkDataSetSurfaceFilter(vtkDeleteDataArrays(ventricles.vol));

    %  check that the surface extracted is actually right - looks fine
    vtkWrite(ventricles.sur,sprintf('%s_ventricles.vtp', cfg));

    % find the RV boundary
    RVFWBoundaryIndices = knnsearch(ventricles.sur.points, RVBoundaryFW.points, 'NSMethod','kdtree');
    RVSeptumBoundaryIndices = knnsearch(ventricles.sur.points, RVBoundarySeptum.points, 'NSMethod','kdtree');
    ventricles.sur.pointData.keep = repmat(uint8(0),size(ventricles.sur.points,1),1);
    indicesToKeepRVFW = unique(RVFWBoundaryIndices);
    indicesToKeepRVSeptum = unique(RVSeptumBoundaryIndices);
    ventricles.sur.pointData.keep(indicesToKeepRVFW) = 1;
    ventricles.sur.pointData.keep(indicesToKeepRVSeptum) = 1;
    RVBoundary = vtkDeleteDataArrays(vtkThreshold(ventricles.sur, 'points', 'keep', [1 1]));
    vtkWrite(RVBoundary,sprintf('%s_RV_boundary.vtp', cfg));

    % find the LV boundary
    LVFWBoundaryIndices = knnsearch(ventricles.sur.points, LVBoundaryFW.points, 'NSMethod','kdtree');
    LVSeptumBoundaryIndices = knnsearch(ventricles.sur.points, LVBoundarySeptum.points, 'NSMethod','kdtree');
    ventricles.sur.pointData.keep = repmat(uint8(0),size(ventricles.sur.points,1),1);
    indicesToKeepLVFW = unique(LVFWBoundaryIndices);
    indicesToKeepLVSeptum = unique(LVSeptumBoundaryIndices);
    ventricles.sur.pointData.keep(indicesToKeepLVFW) = 1;
    ventricles.sur.pointData.keep(indicesToKeepLVSeptum) = 1;
    LVBoundary = vtkDeleteDataArrays(vtkThreshold(ventricles.sur, 'points', 'keep', [1 1]));
    vtkWrite(LVBoundary,sprintf('%s_LV_boundary.vtp', cfg));
end
