function M = baryInterpMat(P, C, targetPoints)
% Computes a matrix to interpolate values from the nodes of a source
% tetrahedral mesh to a set of target points (barycentric interpolation).
%
% M = baryInterpMat(P, C, targetPoints)
%
% Inputs:
%   P: points of the source mesh [numSourcePoints x 3]
%   C: cells of the source mesh [numSourceCells x 4]
%   targetPoints: target points [numTargetPoints x 3]
%
% Outputs:
%   M: interpolation matrix [numTargetPoints x numSourcePoints]
%
% Written by Steffen Schuler, Institute of Biomedical Engineering, KIT

if ~isa(P,'double'), P = double(P); end
if ~isa(C,'double'), C = double(C); end
if ~isa(targetPoints,'double'), targetPoints = double(targetPoints); end

TR = triangulation(C,P);
baryCellIds = pointLocation(TR, targetPoints);

% For points that are not inside any tetrahedron (this may also be points 
% on the boundary of the mesh), find the tetrahedron with the closest
% centroid and use its barycentric coordinates for extrapolation
notFound = isnan(baryCellIds);
if any(notFound)
    centroids = squeeze(mean(reshape(P(C,:),[],size(C,2),size(P,2)),2));
    tmp = vtkMapIds(vtkCreateStruct(centroids), vtkCreateStruct(targetPoints(notFound,:)));
    baryCellIds(notFound) = tmp.pointData.MappedIds;
end

baryPointIds = C(baryCellIds,:);
baryCoords = cartesianToBarycentric(TR, baryCellIds, targetPoints);

m = size(baryCoords,1);
n = size(P,1);
i = reshape(repmat((1:m)',1,4),[],1);
j = baryPointIds(:);
v = baryCoords(:);
M = sparse(i,j,v,m,n);

end