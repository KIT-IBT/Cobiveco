function [vol,mmgOutput] = cobiveco_clipBase(vol, baseNormal, baseOrigin)

baseHeight = baseOrigin(:)'*baseNormal(:);
height = vol.points*baseNormal(:)-baseHeight;

meanEdgLen = mean(vtkEdgeLengths(vol));
mmgSizingParam = [0.1 0.9 1.1];

isovalue = 0;
numTries = 5;
for i = 1:numTries
    [vol,mmgStatus,mmgOutput] = mmg(vol, height, sprintf('-ls %1.5e -nr -hausd %1.5e -hmin %1.5e -hmax %1.5e', isovalue, mmgSizingParam(:)'*meanEdgLen));
    if mmgStatus == 0
        break;
    elseif i < numTries
        warning('Mmg remeshing with an isovalue of %.3f failed (system command status %i). Trying again with a slightly larger isovalue.', isovalue, mmgStatus);
        isovalue = isovalue + 0.1*prctile(max(height(vol.cells),[],2)-min(height(vol.cells),[],2),95);
        mmgSizingParam(1) = 0.8*mmgSizingParam(1);
    else
        warning('Mmg remeshing failed ultimately after trying %i different isovalues (system command status %i).', numTries, mmgStatus);
        return;
    end
end

vol = vtkDeleteDataArrays(vtkThreshold(vol, 'cells', 'class', [2 2]));

end