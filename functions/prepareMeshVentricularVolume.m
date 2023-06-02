function newmesh = prepareMeshVentricularVolume(mesh, original, cfg)

% Prepares the ventricular mesh for tetrahedral volume
%
% newmesh = prepareMeshVentricularVolume(mesh, original, cfg)
%
% Inputs:
%
%   mesh: object of class cobiveco
%   original: object of class cobiveco
%   cfg: configuration, instance in cobiveco class, struct
%
% Output:
% 
%   newmesh: object of class cobiveco (remeshed) 

    disp('Preparing ventricular mesh...');
    t = toc;
    
    isovalue = 0.5;
    numTries = 5;
    for i = 1:numTries
        [newmesh.vol,mmgStatus,mmgOutput2] = mmg(mesh.vol, mesh.ridgeLaplace, sprintf('-ls %1.5e -nr -hausd %1.5e -hmin %1.5e -hmax %1.5e -ar 40', isovalue, cfg.mmgSizingParam(:)'*original.meanEdgLen));
        
        if cfg.exportLevel > 1 || mmgStatus ~= 0
            if i == 1
                fid = fopen(sprintf('%smmgOutput2.txt', cfg.outPrefix), 'w');
            else
                fid = fopen(sprintf('%smmgOutput2.txt', cfg.outPrefix), 'a');
            end
            fprintf(fid, '%s', mmgOutput2);
            fclose(fid);
        end
    
        if mmgStatus == 0
            break;
        elseif i < numTries
            warning('Mmg remeshing with an isovalue of %.3f failed (system command status %i). Trying again with a slightly larger isovalue.', isovalue, mmgStatus);
            isovalue = isovalue + 0.1*prctile(max(mesh.ridgeLaplace(mesh.vol.cells),[],2)-min(mesh.ridgeLaplace(mesh.vol.cells),[],2),95);
            cfg.mmgSizingParam(1) = 0.8*cfg.mmgSizingParam(1);
        else
            error('Mmg remeshing failed ultimately after trying %i different isovalues (system command status %i).', numTries, mmgStatus);
        end
    end
    
    newmesh.sur = vtkDataSetSurfaceFilter(vtkDeleteDataArrays(newmesh.vol));
    newmesh.sur = vtkArrayMapperNearestNeighbor(original.sur, newmesh.sur);
    newmesh.surToVol = vtkMapPointIds(newmesh.vol, newmesh.sur);
    
    P2 = double(newmesh.vol.points);
    C2 = double(newmesh.vol.cells);
    newmesh.L = cotmatrix(P2, C2);
    newmesh.G = grad(P2, C2);
    newmesh.M = baryInterpMat(P2, C2, original.vol.points);
    
    if cfg.exportLevel > 1
        newmesh.debug = newmesh.vol;
        newmesh.debug.pointData.surClass = repmat(uint8(0),size(newmesh.debug.points,1),1);
        newmesh.debug.pointData.surClass(newmesh.surToVol) = newmesh.sur.pointData.class;
    end
    
    
    end