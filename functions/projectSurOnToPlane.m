function [projectedPoints, referenceVectorArray] = projectSurOnToPlane(sur, referenceNormal)
% Project points onto a given plane using a reference reference normal.
% Return in addition an array of vectors caluclate between the projected
% center point and each projected point.
%
% [projectedPoints, referenceVectorArray] = projectSurOnToPlane(sur, referenceNormal)
%
% Inputs:
%
%   sur struct, surface representing the valve annuli with points being [numSurPoints x 3]
%   referenceNormal, double [1x3]
%
% Output:
%
%   projectedPoints, [numSurPoints x 3]
%   referenceVectorArray, double [1x3]
%
% Written by Lisa Pankewitz


% define center of input surface

centerSur = vtkCenterOfArea(sur);

% project points onto the plane of the valve of interest
% center: reference center point,
% determined by projecting surface points onto the valve annulus plane of
% interest
d = (centerSur - sur.points) * referenceNormal';
projectedPoints = sur.points + d .* referenceNormal;

% calculate set of vectors between center and each point in surface

referenceVectorArray = zeros(size(sur.points,1),3);

for n = 1:size(sur.points,1)

    referenceVectorArray(n,:) = (projectedPoints(n,:) - centerSur);
end


end
