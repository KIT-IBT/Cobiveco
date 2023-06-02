function [struct1_OT_boundary, struct1_IT_boundary] = createBoundarySurfacesRT(struct1, OT_surface, IT_surface, cfg)
    % Create 2 boundary surface for the apicobasal coordinate that you can apply readily to the bridge.
    % checks which boundary belong to which site of the bridge (free wall or septum), by comparing it to the reference vector.
    %
    % [struct1_OT_boundary, struct1_IT_boundary] = createBoundarySurfacesRT(struct1, OT_surface, IT_surface, cfg)
    %
    % Inputs:
    %   struct1, struct, tetrahedral vol mesh, with [numVolPoints x 3]
    %   OT_surface, Outflowtract, struct, triangular surface mesh, with [Points x 3]
    %   IT_surface, Inflowtract, struct, triangular surface mesh, with [Points x 3]
    %   cfg configuration prefix for saving the files produced
    %
    % Output:
    %   
    %   struct1_OT_boundary, struct, triangular sur mesh, with [Points x 3]
    %   struct1_IT_boundary, struct, triangular sur mesh, with [Points x 3]
    %
    % Written by Lisa Pankewitz

    % extract surfaces
    struct1.sur= vtkDataSetSurfaceFilter(vtkDeleteDataArrays(struct1.vol));

    %  check that the surface extracted is actually right - looks fine
    vtkWrite(struct1.sur,sprintf('%sstruct1.vtp', cfg));

    % find OT boundary surfaces on struct 1
    struct1_OT_boundary_indices = knnsearch(struct1.sur.points, OT_surface.points, 'NSMethod','kdtree');
    struct1.sur.pointData.OT_boundary = repmat(uint8(0),size(struct1.sur.points,1),1);
    indices_struct1_to_keep_OT = unique(struct1_OT_boundary_indices);
    struct1.sur.pointData.OT_boundary(indices_struct1_to_keep_OT) = 1;
    struct1_OT_boundary  = vtkDeleteDataArrays(vtkThreshold(struct1.sur, 'points', 'OT_boundary', [1 1]));
    vtkWrite(struct1_OT_boundary,sprintf('%sstruct1_OT_boundary.vtp', cfg));

    % find IT boundary surfaces on struct 1
    struct1_IT_boundary_indices = knnsearch(struct1.sur.points, IT_surface.points, 'NSMethod','kdtree');
    struct1.sur.pointData.IT_boundary = repmat(uint8(0),size(struct1.sur.points,1),1);
    indices_struct1_to_keep_IT = unique(struct1_IT_boundary_indices);
    struct1.sur.pointData.IT_boundary(indices_struct1_to_keep_IT) = 1;
    struct1_IT_boundary  = vtkDeleteDataArrays(vtkThreshold(struct1.sur, 'points', 'IT_boundary', [1 1]));
    vtkWrite(struct1_IT_boundary,sprintf('%sstruct1_IT_boundary.vtp', cfg));

end