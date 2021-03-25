function [vol,s,mmgOutput] = cobiveco_removeBridges(vol, sur, baseNormal, ventricles, remesh, thresh, numSubdiv, lambda)

if nargin < 4 || isempty(ventricles)
    ventricles = 'both';
end
if nargin < 5 || isempty(remesh)
    remesh = true;
end
if nargin < 6 || isempty(thresh)
    thresh = 0.1;
end
if nargin < 7 || isempty(numSubdiv)
    numSubdiv = 2;
end
if nargin < 8 || isempty(lambda)
    lambda = 1e2;
end
mmgSizingParam = [0.1 0.9 1.1];

%%
s = vtkLoopSubdivisionFilter(sur, numSubdiv);

if strcmp(ventricles, 'lv')
    endo = vtkThreshold(s, 'points', 'class', [3 3]);
    idsZero = find(s.pointData.class~=3);
elseif strcmp(ventricles, 'rv')
    endo = vtkThreshold(s, 'points', 'class', [4 4]);
    idsZero = find(s.pointData.class~=4);
else % both
    endo = vtkThreshold(s, 'points', 'class', [3 4]);
    idsZero = find(s.pointData.class~=3 & s.pointData.class~=4);
end

baseNormal = baseNormal(:)';
tri = vtkToTriangulation(s);
normals = tri.vertexNormal;
d = normals*baseNormal';
d(idsZero) = 0;
d = 2*max(d-0.5, 0);

edges = vtkFeatureEdges(endo, 1, 0, 0, 0, 0);
idsEdges = vtkMapPointIds(s, edges);

P = double(s.points);
C = double(s.cells);
L = massmatrix(P,C) \ cotmatrix(P,C);
LL = (L'*L);
I = speye(size(L));

d = (I+lambda*LL) \ d;
d(idsZero) = 0;
d(idsEdges) = 100*d(idsEdges);

d = (I+10*lambda*LL) \ d;
d(idsZero) = 0;

s.pointData.d = d;
bridgeSur = vtkThreshold(s, 'points', 'd', [thresh inf]);
bridgeSur = vtkDeleteDataArrays(bridgeSur);

el = mean(vtkEdgeLengths(bridgeSur));
n = 100;
bridgeGrid = repmat(bridgeSur.points, n, 1) - repmat(reshape(repmat((-1:n-2)*el, size(bridgeSur.points,1), 1), [], 1), 1, 3) .* repmat(baseNormal,n*size(bridgeSur.points,1),1);

%%
idsBridgeGrid = knnsearch(bridgeGrid, vol.points);

dist = sqrt(sum((bridgeGrid(idsBridgeGrid,:)-vol.points).^2,2));
subId = find(dist < el);
projPoints = vol.points(subId,:)-(vol.points(subId,:)*baseNormal')*baseNormal;

bridgeSurFine = vtkLinearSubdivisionFilter(bridgeSur,1);
bridgeSurFine.points = bridgeSurFine.points-(bridgeSurFine.points*baseNormal')*baseNormal;

idsBridgeSurFine = knnsearch(bridgeSurFine.points, projPoints);
dist(subId) = sqrt(sum((bridgeSurFine.points(idsBridgeSurFine,:)-projPoints).^2,2));

vol.pointData.dist = dist;

%%
if remesh
    meanEdgLen = mean(vtkEdgeLengths(vol));
    isovalue = el;
    numTries = 5;
    for i = 1:numTries
        [vol,mmgStatus,mmgOutput] = mmg(vol, dist, sprintf('-ls %1.5e -nr -hausd %1.5e -hmin %1.5e -hmax %1.5e', el, mmgSizingParam(:)'*meanEdgLen));
        if mmgStatus == 0
            break;
        elseif i < numTries
            warning('Mmg remeshing with an isovalue of %.3f failed (system command status %i). Trying again with a slightly larger isovalue.', isovalue, mmgStatus);
            isovalue = isovalue + 0.1*el;
            mmgSizingParam(1) = 0.8*mmgSizingParam(1);
        else
            warning('Mmg remeshing failed ultimately after trying %i different isovalues (system command status %i).', numTries, mmgStatus);
            return;
        end
    end
    vol = vtkThreshold(vol, 'cells', 'class', [2 2]);
    vol = vtkConnectivityFilter(vtkDeleteDataArrays(vol));
    rid = double(mode(vol.cellData.RegionId));
    vol = vtkDeleteDataArrays(vtkThreshold(vol, 'cells', 'RegionId', [rid rid]));
else
    mmgOutput = '';
end

end