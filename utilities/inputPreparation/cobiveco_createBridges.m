function [vol,mmgOutput] = cobiveco_clipBridges(mesh, baseNormal_rv, baseOrigin_rv, baseNormal_lv, baseOrigin_lv)
    % script to create bridges for each ventricle
    % mesh needs to have the tv information as pointData
    % we will clip the rv and lv bridges individually to just get rid of the bridges
    % defining the clipping height for the LV
    baseHeight_lv = baseOrigin_lv(:)'*baseNormal_lv(:);
    % defining the clipping height for the RV
    baseHeight_rv = baseOrigin_rv(:)'*baseNormal_rv(:);
    
    % Cutting  plane (not curved)
    height = zeros(size(mesh.vol.points,1),1);

    for i = 1:size(mesh.tv)
        % if we are in the rv, we apply the cutoff for the rv bridge including the tricuspid and pulmonary valve annuli
        if mesh.tv(i) == 1
            height(i) = mesh.vol.points(i,:)*baseNormal_rv(:)-baseHeight_rv;
        else
            %else mesh.tv(i) == 0, meaning we are in the lv
            height(i) = mesh.vol.points(i,:)*baseNormal_lv(:)-baseHeight_lv;
        end
    end
    

    meanEdgLen = mean(vtkEdgeLengths(mesh.vol));
    mmgSizingParam = [0.1 0.9 1.1];
    
    isovalue = 0;
    numTries = 5;
    for i = 1:numTries
        [vol,mmgStatus,mmgOutput] = mmg(mesh.vol, height, sprintf('-ls %1.5e -nr -hausd %1.5e -hmin %1.5e -hmax %1.5e', isovalue, mmgSizingParam(:)'*meanEdgLen));
        %keyboard
        if mmgStatus == 0
            break;
        elseif i < numTries
            warning('Mmg remeshing with an isovalue of %.3f failed (system command status %i). Trying again with a slightly larger isovalue.', isovalue, mmgStatus);
            %keyboard
            isovalue = isovalue + 0.1*prctile(max(height(mesh.vol.cells),[],2)-min(height(mesh.vol.cells),[],2),95);
            mmgSizingParam(1) = 0.8*mmgSizingParam(1);
        else
            warning('Mmg remeshing failed ultimately after trying %i different isovalues (system command status %i).', numTries, mmgStatus);
            return;
        end
    end
    
    % create bridges for lv and rv
    vol = vtkDeleteDataArrays(vtkThreshold(vol, 'cells', 'class', [2 2]));
    
    end