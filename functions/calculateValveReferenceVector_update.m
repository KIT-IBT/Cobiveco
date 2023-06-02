function [referenceVectorMatch, referenceVectorProjectedPlane, referenceVectorArray, idxPointOnSur] = calculateValveReferenceVectorUpdate(sur, referenceNormal, referenceCenter)
% Calculate the vectors connecting the center of a surface
% to each point in the surface and find the one vector matching 
% a reference vector defined best by minimzing the angle.
%
% [referenceVectorMatch, referenceVectorProjectedPlane, referenceVectorArray] = calculateValveReferenceVectorUpdate(sur, referenceNormal, referenceCenter)
%
% Inputs:
%
%   sur, surface structure representing the valve annuli with points being [numSurPoints x 3]
%   referenceVector, double [1x3]
%   normal, double [1x3]
%
% Output:
%
%   referenceVectorMatch, vector matching the referencevectors best by minimizing the angle, double [1x3]
%
%
% Written by Lisa Pankewitz

% define center of input surface

centerSur = vtkCenterOfArea(sur);

% project points onto the plane of the valve of interest
% center: reference center point,
% determined by projecting surface points onto the valve annulus plane of
% interest
d = (referenceCenter - sur.points) * referenceNormal';
projectedPoints = sur.points + d .* referenceNormal;

% project center, too
dCenter = (referenceCenter - centerSur) * referenceNormal';
projectedCenter = centerSur + dCenter .* referenceNormal;

% calculate set of vectors between center and each point in surface
referenceVectorArray = zeros(size(sur.points,1),3);

for n = 1:size(sur.points,1)

    referenceVectorArray(n,:) = (projectedPoints(n,:) - projectedCenter);
end

% calculate reference vector refVectorToReferenceCenter
refVectorToReferenceCenter = referenceCenter -projectedCenter;

% find the one that aligns best with the refVectorToReferenceCenter
% compare each vector to the referencecevtor
angleVectorToMvCenter = zeros(size(sur.points,1),1);

for n = 1:size(sur.points,1)
    angleVectorToReferenceVector(n) = atan2(norm(cross(referenceVectorArray(n,:)',refVectorToReferenceCenter')),dot(referenceVectorArray(n,:)',refVectorToReferenceCenter'));
end


% find the vector which produces the smallest angle compared to the referencevector
[value,idxPointOnSur] = min(abs(angleVectorToReferenceVector));

% calculate set of vectors between center and each point in surface
referenceVectorArrayNotProjected = zeros(size(sur.points,1),3);

for n = 1:size(sur.points,1)

    referenceVectorArrayNotProjected(n,:) = (sur.points(n,:) - centerSur);
end


referenceVectorMatch = referenceVectorArrayNotProjected(idxPointOnSur,:);

referenceVectorProjectedPlane = referenceVectorArray(idxPointOnSur,:);

end
