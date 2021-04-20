function prepareMesh2(o)

if ~o.available.transmural
    o.computeTransmural;
end
o.printStatus('Preparing mesh 2...');
t = toc;

isovalue = 0.5;
numTries = 5;
for i = 1:numTries
    [o.m2.vol,mmgStatus,mmgOutput2] = mmg(o.m1.vol, o.m1.ridgeLaplace, sprintf('-ls %1.5e -nr -hausd %1.5e -hmin %1.5e -hmax %1.5e', isovalue, o.cfg.mmgSizingParam(:)'*o.m0.meanEdgLen));

    if o.cfg.exportLevel > 1 || mmgStatus ~= 0
        if i == 1
            fid = fopen(sprintf('%smmgOutput2.txt', o.cfg.outPrefix), 'w');
        else
            fid = fopen(sprintf('%smmgOutput2.txt', o.cfg.outPrefix), 'a');
        end
        fprintf(fid, '%s', mmgOutput2);
        fclose(fid);
    end

    if mmgStatus == 0
        break;
    elseif i < numTries
        warning('Mmg remeshing with an isovalue of %.3f failed (system command status %i). Trying again with a slightly larger isovalue.', isovalue, mmgStatus);
        isovalue = isovalue + 0.1*prctile(max(o.m1.ridgeLaplace(o.m1.vol.cells),[],2)-min(o.m1.ridgeLaplace(o.m1.vol.cells),[],2),95);
        o.cfg.mmgSizingParam(1) = 0.8*o.cfg.mmgSizingParam(1);
    else
        error('Mmg remeshing failed ultimately after trying %i different isovalues (system command status %i).', numTries, mmgStatus);
    end
end

o.m2.sur = vtkDataSetSurfaceFilter(vtkDeleteDataArrays(o.m2.vol));
o.m2.sur = vtkArrayMapperNearestNeighbor(o.m0.sur, o.m2.sur);
o.m2.surToVol = vtkMapPointIds(o.m2.vol, o.m2.sur);

P2 = double(o.m2.vol.points);
C2 = double(o.m2.vol.cells);
o.m2.L = cotmatrix(P2, C2);
o.m2.G = grad(P2, C2);
o.m2.M = baryInterpMat(P2, C2, o.m0.vol.points);

if o.cfg.exportLevel > 1
    o.m2.debug = o.m2.vol;
    o.m2.debug.pointData.surClass = repmat(uint8(0),size(o.m2.debug.points,1),1);
    o.m2.debug.pointData.surClass(o.m2.surToVol) = o.m2.sur.pointData.class;
end

o.printStatus(sprintf('%.1f seconds\n', toc-t), true);
o.available.mesh2 = true;

end