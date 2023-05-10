function [ventrvol,rv_bridge_vol, lv_bridge_vol, mmgOutput1,  mmgOutput2,  mmgOutput3, idx_rv_bridge,idx_lv_bridge,ventricle_filter] = cobiveco_clipBridges(mesh, idx_lv_bridge, idx_rv_bridge)
    %
    % Remeshing the intravalvular regions based on a filter given by node indices from the original mesh.
    % Output contains three tetrahedral volume meshes.
    %
    % [ventrvol,rv_bridge_vol, lv_bridge_vol, mmgOutput1,  mmgOutput2,  mmgOutput3, idx_rv_bridge,idx_lv_bridge,ventricle_filter] = cobiveco_clipBridges(mesh, idx_lv_bridge, idx_rv_bridge)
    % Inputs:
    %   mesh, struct, tetrahedral vol mesh, with [numVolPoints x 3]
    %   idx_lv_bridge , int [nx3], with n being number of nodes in the LV bridge
    %   idx_rv_bridge , int [nx3], with n being number of nodes in the RV bridge
    %
    % Output:
    %
    %   ventrvol, output struct, a tetrahedral mesh, with [numVolPointsRVBridge x 3]
    %   rv_bridge_vol, output struct, a tetrahedral mesh, with [numVolPointsRVBridge x 3]
    %   rv_bridge_vol, output struct, a tetrahedral mesh, with [numVolPointsRVBridge x 3]
    %
    % Written by Lisa Pankewitz

    %rv_bridge_filter = ones(size(mesh.vol.points,1),1);
    %rv_bridge_filter = -0.001*rv_bridge_filter;
    %rv_bridge_filter(idx_rv_bridge) = 1;
    rv_bridge_filter = zeros(size(mesh.vol.points,1),1);
    rv_bridge_filter(idx_rv_bridge) = 1;

    meanEdgLen = mean(vtkEdgeLengths(mesh.vol));
    mmgSizingParam = [0.1 0.9 1.1];

    isovalue = 0;
    numTries = 5;
    for i = 1:numTries
        [vol,mmgStatus,mmgOutput2] = mmg(mesh.vol, rv_bridge_filter, sprintf('-ls %1.5e -nr -hausd %1.5e -hmin %1.5e -hmax %1.5e', isovalue, mmgSizingParam(:)'*meanEdgLen));
        %keyboard
        if mmgStatus == 0
            break;
        elseif i < numTries
            warning('Mmg remeshing with an isovalue of %.3f failed (system command status %i). Trying again with a slightly larger isovalue.', isovalue, mmgStatus);
            %keyboard
            isovalue = isovalue + 0.1*prctile(max(rv_bridge_filter(mesh.vol.cells),[],2)-min(rv_bridge_filter(mesh.vol.cells),[],2),95);
            mmgSizingParam(1) = 0.8*mmgSizingParam(1);
        else
            warning('Mmg remeshing failed ultimately after trying %i different isovalues (system command status %i).', numTries, mmgStatus);
            return;
        end
    end

    % create rv bridge
    rv_bridge_vol = vtkDeleteDataArrays(vtkThreshold(vol, 'cells', 'class', [2 2]));

    disp('Finished rv bridge.');

    % alternative code no for loops
    ventricle_filter = ones(size(mesh.vol.points,1),1);
    ventricle_filter(idx_rv_bridge) = -0.01;
    ventricle_filter(idx_lv_bridge) = -0.01;

    %returns elapsed time
    t = toc;
    meanEdgLen = mean(vtkEdgeLengths(mesh.vol));
    mmgSizingParam = [0.1 0.9 1.1];

    isovalue = 0;
    numTries = 5;
    for i = 1:numTries
        [vol,mmgStatus,mmgOutput1] = mmg(mesh.vol, ventricle_filter, sprintf('-ls %1.5e -nr -hausd %1.5e -hmin %1.5e -hmax %1.5e', isovalue, mmgSizingParam(:)'*meanEdgLen));

        if mmgStatus == 0
            break;
        elseif i < numTries
            warning('Mmg remeshing with an isovalue of %.3f failed (system command status %i). Trying again with a slightly larger isovalue.', isovalue, mmgStatus);

            isovalue = isovalue + 0.1*prctile(max(ventricle_filter(mesh.vol.cells),[],2)-min(ventricle_filter(mesh.vol.cells),[],2),95);
            mmgSizingParam(1) = 0.8*mmgSizingParam(1);
        else
            warning('Mmg remeshing failed ultimately after trying %i different isovalues (system command status %i).', numTries, mmgStatus);
            return;
        end
    end

    % returns elapsed time
    t = toc;
    % create ventricular volume
    ventrvol = vtkDeleteDataArrays(vtkThreshold(vol, 'cells', 'class', [2 2]));
    disp('Ventricular Mesh done')

    % filter only lv bridge
    lv_bridge_filter = zeros(size(mesh.vol.points,1),1);
    lv_bridge_filter(idx_lv_bridge) = 1;


    meanEdgLen = mean(vtkEdgeLengths(mesh.vol));
    mmgSizingParam = [0.1 0.9 1.1];

    isovalue = 0;
    numTries = 5;
    for i = 1:numTries
        [vol,mmgStatus,mmgOutput3] = mmg(mesh.vol, lv_bridge_filter, sprintf('-ls %1.5e -nr -hausd %1.5e -hmin %1.5e -hmax %1.5e', isovalue, mmgSizingParam(:)'*meanEdgLen));
        %keyboard
        if mmgStatus == 0
            break;
        elseif i < numTries
            warning('Mmg remeshing with an isovalue of %.3f failed (system command status %i). Trying again with a slightly larger isovalue.', isovalue, mmgStatus);
            %keyboard
            isovalue = isovalue + 0.1*prctile(max(lv_bridge_filter(mesh.vol.cells),[],2)-min(lv_bridge_filter(mesh.vol.cells),[],2),95);
            mmgSizingParam(1) = 0.8*mmgSizingParam(1);
        else
            warning('Mmg remeshing failed ultimately after trying %i different isovalues (system command status %i).', numTries, mmgStatus);
            return;
        end
    end

    disp('LV bridge done')
    % create lv bridge

    lv_bridge_vol = vtkDeleteDataArrays(vtkThreshold(vol, 'cells', 'class', [2 2]));

    end
