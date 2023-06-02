function [sur,debug] = cobiveco_createClasses(sur, baseNormal, maxAngle, numSubdiv)

if nargin < 4 || isempty(numSubdiv)
    numSubdiv = 0;
end

% figure out points that belong to rv or lv
ids_lv = find(sur.tv == 1);
ids_rv = find(sur.tv == 0);
% remove rv or lv from struct
sur.pointData.tv = sur.tv;
surLv = vtkDataSetSurfaceFilter(vtkThreshold(sur, 'points', 'tv', [1 1]));
surRv = vtkDataSetSurfaceFilter(vtkThreshold(sur, 'points', 'tv', [0 0]));
s_Lv = vtkLoopSubdivisionFilter(surLv, numSubdiv);
tr_Lv = vtkToTriangulation(s_Lv);
dotProd_Lv = double(tr_Lv.faceNormal * baseNormal');
dotProd_Lv = min(max(dotProd_Lv,-1),1);
s_Lv.cellData.angle = acos(-dotProd_Lv)/pi*180;

% find centroids 
centroids_Lv = vtkCellCentroids(s_Lv);
center_Lv = vtkCenterOfArea(s_Lv);
s_Lv.cellData.centerDist = (centroids_Lv.points-center_Lv)*baseNormal';
s_Lv.pointData.centerDist = (s_Lv.points-center_Lv)*baseNormal';

s_Lv.cellData.base = double(s_Lv.cellData.angle < maxAngle & s_Lv.cellData.centerDist < 0);
s_Lv.pointData.ids = int32(1:size(s_Lv.points,1))';

base_Lv = vtkConnectivityFilter(vtkThreshold(s_Lv, 'cells', 'base', [1 inf]));
idsUnknown = [];
if numel(unique(base_Lv.pointData.RegionId)) > 1
    rid = mode(base_Lv.pointData.RegionId);
    idsUnknown = base_Lv.pointData.ids(base_Lv.pointData.RegionId ~= rid);
end

nonbase_Lv = vtkConnectivityFilter(vtkThreshold(s_Lv, 'cells', 'base', [-inf 0]));
if numel(unique(nonbase_Lv.pointData.RegionId)) ~= 3
    warning('cobiveco_createClasses: Number of non-base regions ~= 3. Manual interaction needed - aborting.');
    debug = s_Lv;
    return;
end
region = cell(3,1);
radius = NaN(3,1);
region{1} = vtkThreshold(nonbase_Lv, 'points', 'RegionId', [0 0]);
region{2} = vtkThreshold(nonbase_Lv, 'points', 'RegionId', [1 1]);
region{3} = vtkThreshold(nonbase_Lv, 'points', 'RegionId', [2 2]);
[~,radius(1)] = vtkSmallestEnclosingSphere(region{1});
[~,radius(2)] = vtkSmallestEnclosingSphere(region{2});
[~,radius(3)] = vtkSmallestEnclosingSphere(region{3});
[~,epiInd] = max(radius);
endoInd = setdiff(1:3,epiInd);
epi = region{epiInd};
endo1 = region{endoInd(1)};
endo2 = region{endoInd(2)};
contour1 = vtkContourFilter(endo1, 'points', 'centerDist', 0);
contour2 = vtkContourFilter(endo2, 'points', 'centerDist', 0);
d1 = min(sqrt(sum((contour1.points-center_Lv).^2,2)));
d2 = min(sqrt(sum((contour2.points-center_Lv).^2,2)));

ids_epi = vtkMapPointIds(s_Lv, epi);
if d1 < d2
    ids_lv = vtkMapPointIds(s_Lv, endo1);
    ids_rv = vtkMapPointIds(s_Lv, endo2);
else
    ids_lv = vtkMapPointIds(s_Lv, endo2);
    ids_rv = vtkMapPointIds(s_Lv, endo1);
end

s_Lv.pointData.class = ones(size(s_Lv.points,1), 1, 'uint8');
s_Lv.pointData.class(ids_epi) = 2;
s_Lv.pointData.class(ids_lv) = 3;
s_Lv.pointData.class(ids_rv) = 4;

if any(idsUnknown)
    L = cotmatrix(double(s_Lv.points), double(s_Lv.cells));
    idsKnown = setdiff(1:size(s_Lv.points,1), idsUnknown);
    s_Lv.pointData.class = uint8(round(solveLaplace(L, idsKnown, double(s_Lv.pointData.class(idsKnown)))));
end

debug = s_Lv;
s_Lv = vtkDeleteDataArrays(s_Lv, {'pointData.class'});

if numSubdiv > 0
    sur.pointData.class = zeros(size(sur.points,1), 1, 'uint8');
    ids = vtkMapPointIds(sur, s_Lv);
    for i = 1:size(sur.points,1)
        sur.pointData.class(i) = mode(s_Lv.pointData.class(ismember(ids, i)));
    end
else
    sur = s_Lv;
end

vtkWrite(sur,'result_gmsh_split20211119/testsur.vtp');
end