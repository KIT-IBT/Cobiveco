function [ventrvol,rv_bridge_vol, lv_bridge_vol, idx_rv_bridge, idx_lv_bridge, ventricle_filter] = cobiveco_clipBridges_noremesh(mesh, baseNormal_rv, baseOrigin_rv, filter_lv_bridge, original)
    % mesh needs to have the tv information as pointData
    % we will clip the rv and lv bridges individually to just get rid of the bridges
    %we will just use vtk delete array to get rid of the bridges
    % defining the clipping height for the RV
    baseHeight_rv = baseOrigin_rv(:)'*baseNormal_rv(:);
    
    % Cutting  plane (not curved)
    height = zeros(size(mesh.vol.points,1),1);

    %returns elapsed time
    t = toc;

    meanEdgLen = original.meanEdgLen;

    % create bridges
    % start with rv bridge
    % Cutting  plane (not curved)

    baseHeight_rv = baseOrigin_rv(:)'*-baseNormal_rv(:);
    

    height_rv = zeros(size(mesh.vol.points,1),1);
    
    for i = 1:size(mesh.tv)
        % if we are in the rv, we apply the cutoff for the rv bridge including the tricuspid and pulmonary valve annuli
        if mesh.tv(i) == 1
            height_rv(i) = mesh.vol.points(i,:)*(-baseNormal_rv(:))-baseHeight_rv - height(i);
        else
            %else mesh.tv(i) == 0, meaning we are in the lv
            height_rv(i) = -100;
        end
    end

    idx_rv_bridge = find(height_rv > 0);
    %ventricles 0
    ventricle_filter = zeros(size(mesh.vol.points,1),1);
    % rv 1
    ventricle_filter(idx_rv_bridge) = 1;
    idx_lv_bridge = find(filter_lv_bridge == 1);
    %lv 2
    ventricle_filter(idx_lv_bridge) = 2;


    mesh.vol.pointData.bridgefilter = repmat(uint8(0),size(mesh.vol.points,1),1);
    mesh.vol.pointData.bridgefilter = ventricle_filter;
    % need to figure out the elements and then based on the elements back calculate filter 

    % create rv bridge
    rv_bridge_vol = vtkDeleteDataArrays(vtkThreshold(mesh.vol, 'points', 'bridgefilter', [1 1]));

    % create ventricular volume
    ventrvol = vtkDeleteDataArrays(vtkThreshold(mesh.vol, 'points', 'bridgefilter', [0 0]));

    % create lv bridge
    lv_bridge_vol = vtkDeleteDataArrays(vtkThreshold(mesh.vol, 'points', 'bridgefilter', [2 2]));

    end