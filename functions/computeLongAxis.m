function longAx = computeLongAxis(sur, directionVec)
% Computes the heart's long axis as the vector minimizing the dot product
% with the surface normals of the LV endocardium, i.e. the vector that is
% "most orthogonal" to the surface normals.
%
% longAx = computeLongAxis(sur, basePoints)
%
% Inputs:
%   sur: LV endocardial surface as VTK struct
%   directionVec: A vector coarsely directed from base towards apex
%
% Outputs:
%   longAx: Unit vector directed from base towards apex
%
% Written by Steffen Schuler, Institute of Biomedical Engineering, KIT

TR = vtkToTriangulation(sur);
normals = TR.faceNormal;
areaWeights = doublearea(TR.Points, TR.ConnectivityList);
areaWeights = areaWeights/mean(areaWeights);

% p-norm is chosen to have average properties of 1-norm and 2-norm
h = (1/sqrt(2)+1)/2;
p = 1/(log2(sqrt(2)/h));

objFun = @(longAxis) sum(abs(areaWeights.*normals*longAxis').^p) + size(normals,1)*abs(norm(longAxis)-1)^p;
options = optimset('MaxFunEvals', 1e4);
longAx = NaN(3);
objVal = NaN(3,1);
[longAx(1,:),objVal(1)] = fminsearch(objFun, [1 0 0], options);
[longAx(2,:),objVal(2)] = fminsearch(objFun, [0 1 0], options);
[longAx(3,:),objVal(3)] = fminsearch(objFun, [0 0 1], options);
[~,minInd] = min(objVal);
longAx = longAx(minInd,:);
longAx = longAx/norm(longAx);

% make sure longAxis is directed from base towards apex
d = directionVec(:)' * longAx';
longAx = sign(d) * longAx;

end