function prepareMesh1(o)

if ~o.available.transventricular
    o.computeTransventricular;
end
o.printStatus('Preparing mesh 1...');
t = toc;

isovalue = 0.5;
numTries = 5;
for i = 1:numTries
    [o.m1.vol,mmgStatus,mmgOutput1] = mmg(o.m0.vol, o.m0.tvLaplace, sprintf('-ls %1.5e -nr -hausd %1.5e -hmin %1.5e -hmax %1.5e', isovalue, o.cfg.mmgSizingParam(:)'*o.m0.meanEdgLen));

    if o.cfg.exportLevel > 1 || mmgStatus ~= 0
        if i == 1
            fid = fopen(sprintf('%smmgOutput1.txt', o.cfg.outPrefix), 'w');
        else
            fid = fopen(sprintf('%smmgOutput1.txt', o.cfg.outPrefix), 'a');
        end
        fprintf(fid, '%s', mmgOutput1);
        fclose(fid);
    end

    if mmgStatus == 0
        break;
    elseif i < numTries
        warning('Mmg remeshing with an isovalue of %.3f failed (system command status %i). Trying again with a slightly larger isovalue.', isovalue, mmgStatus);
        isovalue = isovalue + 0.1*prctile(max(o.m0.tvLaplace(o.m0.vol.cells),[],2)-min(o.m0.tvLaplace(o.m0.vol.cells),[],2),95);
        o.cfg.mmgSizingParam(1) = 0.8*o.cfg.mmgSizingParam(1);
    else
        error('Mmg remeshing failed ultimately after trying %i different isovalues (system command status %i).', numTries, mmgStatus);
    end
end

o.m1.sur = vtkDataSetSurfaceFilter(vtkDeleteDataArrays(o.m1.vol));
o.m1.sur = vtkArrayMapperNearestNeighbor(o.m0.sur, o.m1.sur);
o.m1.surToVol = vtkMapPointIds(o.m1.vol, o.m1.sur);
o.m1.meanEdgLen = mean(vtkEdgeLengths(o.m1.vol));

P1 = double(o.m1.vol.points);
C1 = double(o.m1.vol.cells);
o.m1.L = cotmatrix(P1, C1);
o.m1.G = grad(P1, C1);
o.m1.massMat = massmatrix(P1, C1, 'voronoi');
o.m1.M = baryInterpMat(P1, C1, o.m0.vol.points);

if o.cfg.exportLevel > 1
    o.m1.debug = o.m1.vol;
    o.m1.debug.pointData.surClass = repmat(uint8(0),size(o.m1.debug.points,1),1);
    o.m1.debug.pointData.surClass(o.m1.surToVol) = o.m1.sur.pointData.class;
end

o.printStatus(sprintf('%.1f seconds\n', toc-t), true);
o.available.mesh1 = true;

end