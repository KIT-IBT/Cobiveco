function [ventricularmesh, rv_bridge, lv_bridge] = SplitMesh(mesh, cfg)
    %
    % [idx_non_bridge] = defineBridgesShortestPath(mesh,path1, path2, epi_path1, epi_path2, tv)
    % Inputs:
    %   mesh, struct, tehtrahedral vol mesh, with [numVolPoints x 3] and field RegionIdVentricle
    %   cfg configuration struct with output prefix
    %
    % Output:
    %   
    %   ventricularmesh, output struct, a tetrahedral mesh, with [numVolPointsRVBridge x 3]
    %   rv_bridge, output struct, a tetrahedral mesh, with [numVolPointsRVBridge x 3]
    %   lv_bridge, output struct, a tetrahedral mesh, with [numVolPointsRVBridge x 3]
    %
    % Written by Lisa Pankewitz

    try
        isfield(mesh.vol.cellData,'RegionIdVentricle') && ~isempty(mesh.vol.cellData.RegionIdVentricle)

    catch
        error('Field RegionIdVentricle does not exist.\n')
    end

    lv_bridge.vol = vtkDeleteDataArrays(vtkThreshold(mesh.vol, 'cells', 'RegionIdVentricle', [2 2]));
    rv_bridge.vol = vtkDeleteDataArrays(vtkThreshold(mesh.vol, 'cells', 'RegionIdVentricle', [1 1]));
    ventricularmesh.vol = vtkDeleteDataArrays(vtkThreshold(mesh.vol, 'cells', 'RegionIdVentricle', [0 0]));
    % cell2pointdata
    outStruct = vtkCellDataToPointData(mesh.vol, (true));

    % check
    vtkWrite(rv_bridge.vol, sprintf('%srv_bridge.vtu', cfg.outPrefix));
    vtkWrite(lv_bridge.vol, sprintf('%slv_bridge.vtu', cfg.outPrefix));
    vtkWrite(ventricularmesh.vol, sprintf('%sventricularmesh.vtu', cfg.outPrefix));

end