function newmesh = prepareSeptalSurface(mesh, cfg)

% Remeshing step to allow for extraction of septal surface.
%
% newmesh = prepareSeptalSurface(mesh, cfg)
%
% Inputs:
%
%   mesh: object of class cobiveco, tetrahedral mesh
%   cfg: configuration instance in cobiveco class, struct
%
% Output:
%
%   newmesh: object of class cobiveco (remeshed)

disp('Preparing septal surface...');
t = toc;

isovalue = 0.5;
numTries = 5;
for i = 1:numTries
    [newmesh.vol,mmgStatus,mmgOutput1] = mmg(mesh.vol, mesh.tvLaplace, sprintf('-ls %1.5e -nr -hausd %1.5e -hmin %1.5e -hmax %1.5e', isovalue, cfg.mmgSizingParam(:)'*mesh.meanEdgLen));

    if cfg.exportLevel > 1 || mmgStatus ~= 0
        if i == 1
            fid = fopen(sprintf('%smmgOutput1.txt', cfg.outPrefix), 'w');
        else
            fid = fopen(sprintf('%smmgOutput1.txt', cfg.outPrefix), 'a');
        end
        fprintf(fid, '%s', mmgOutput1);
        fclose(fid);
    end

    if mmgStatus == 0
        break;
    elseif i < numTries
        warning('Mmg remeshing with an isovalue of %.3f failed (system command status %i). Trying again with a slightly larger isovalue.', isovalue, mmgStatus);
        isovalue = isovalue + 0.2*prctile(max(mesh.tvLaplace(mesh.vol.cells),[],2)-min(mesh.tvLaplace(mesh.vol.cells),[],2),95);
        cfg.mmgSizingParam(1) = 1.3*cfg.mmgSizingParam(1);
    else
        error('Mmg remeshing failed ultimately after trying %i different isovalues (system command status %i).', numTries, mmgStatus);
    end
end

newmesh.sur = vtkDataSetSurfaceFilter(vtkDeleteDataArrays(newmesh.vol));
newmesh.sur = vtkArrayMapperNearestNeighbor(mesh.sur, newmesh.sur);
newmesh.surToVol = vtkMapPointIds(newmesh.vol, newmesh.sur);
newmesh.meanEdgLen = mean(vtkEdgeLengths(newmesh.vol));

P1 = double(newmesh.vol.points);
C1 = double(newmesh.vol.cells);
newmesh.L = cotmatrix(P1, C1);
newmesh.G = grad(P1, C1);
newmesh.massMat = massmatrix(P1, C1, 'voronoi');
newmesh.M = baryInterpMat(P1, C1, mesh.vol.points);

if cfg.exportLevel > 1
    newmesh.debug = newmesh.vol;
    newmesh.debug.pointData.surClass = repmat(uint8(0),size(newmesh.debug.points,1),1);
    newmesh.debug.pointData.surClass(newmesh.surToVol) = newmesh.sur.pointData.class;
end



end
