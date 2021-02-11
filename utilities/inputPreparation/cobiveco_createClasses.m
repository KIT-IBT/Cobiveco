function [sur,debug] = cobiveco_createClasses(sur, baseNormal, maxAngle, numSubdiv)

if nargin < 4 || isempty(numSubdiv)
    numSubdiv = 0;
end

s = vtkLoopSubdivisionFilter(sur, numSubdiv);
tr = vtkToTriangulation(s);
dotProd = double(tr.faceNormal * baseNormal');
dotProd = min(max(dotProd,-1),1);
s.cellData.angle = acos(-dotProd)/pi*180;

centroids = vtkCellCentroids(s);
center = vtkCenterOfArea(s);
s.cellData.centerDist = (centroids.points-center)*baseNormal';
s.pointData.centerDist = (s.points-center)*baseNormal';

s.cellData.base = double(s.cellData.angle < maxAngle & s.cellData.centerDist < 0);
s.pointData.ids = int32(1:size(s.points,1))';

base = vtkConnectivityFilter(vtkThreshold(s, 'cells', 'base', [1 inf]));
idsUnknown = [];
if numel(unique(base.pointData.RegionId)) > 1
    rid = mode(base.pointData.RegionId);
    idsUnknown = base.pointData.ids(base.pointData.RegionId ~= rid);
end

nonbase = vtkConnectivityFilter(vtkThreshold(s, 'cells', 'base', [-inf 0]));
if numel(unique(nonbase.pointData.RegionId)) ~= 3
    warning('cobiveco_createClasses: Number of non-base regions ~= 3. Manual interaction needed - aborting.');
    debug = s;
    return;
end
region = cell(3,1);
radius = NaN(3,1);
region{1} = vtkThreshold(nonbase, 'points', 'RegionId', [0 0]);
region{2} = vtkThreshold(nonbase, 'points', 'RegionId', [1 1]);
region{3} = vtkThreshold(nonbase, 'points', 'RegionId', [2 2]);
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
d1 = min(sqrt(sum((contour1.points-center).^2,2)));
d2 = min(sqrt(sum((contour2.points-center).^2,2)));

ids_epi = vtkMapPointIds(s, epi);
if d1 < d2
    ids_lv = vtkMapPointIds(s, endo1);
    ids_rv = vtkMapPointIds(s, endo2);
else
    ids_lv = vtkMapPointIds(s, endo2);
    ids_rv = vtkMapPointIds(s, endo1);
end

s.pointData.class = ones(size(s.points,1), 1, 'uint8');
s.pointData.class(ids_epi) = 2;
s.pointData.class(ids_lv) = 3;
s.pointData.class(ids_rv) = 4;

if any(idsUnknown)
    L = cotmatrix(double(s.points), double(s.cells));
    idsKnown = setdiff(1:size(s.points,1), idsUnknown);
    s.pointData.class = uint8(round(solveLaplace(L, idsKnown, double(s.pointData.class(idsKnown)))));
end

debug = s;
s = vtkDeleteDataArrays(s, {'pointData.class'});

if numSubdiv > 0
    sur.pointData.class = zeros(size(sur.points,1), 1, 'uint8');
    ids = vtkMapPointIds(sur, s);
    for i = 1:size(sur.points,1)
        sur.pointData.class(i) = mode(s.pointData.class(ismember(ids, i)));
    end
else
    sur = s;
end

end