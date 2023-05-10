function referenceVectorMatch = calculateValveReferenceVector(sur, referenceVector)
% Calculate the vectors connecting the center of a surface
% to each point in the surface and find the one vector matching
% a referenceVector defined best by minimzing the angle.
%
% referenceVectorMatch = calculateValveReferenceVector(sur, referenceVector)
%
% Inputs:
%
%   sur struct, surface representing the valve annuli with points being [numSurPoints x 3]
%   referenceVector, double [1x3]
%
% Output:
%
%   referenceVectorMatch, Vector matching the referenceVector best by minimizing the angle, double [1x3]
%
%
% Written by Lisa Pankewitz


% define center of input surface

centerSur = vtkCenterOfArea(sur);

% calculate set of vectors between center and each point in surface

referenceVectorArray = zeros(size(sur.points,1),3);

for n = 1:size(sur.points,1)

    referenceVectorArray(n,:) = (sur.points(n,:) - centerSur);
end

% compare each vector to the reference vector

angleVectorToReferenceVector = zeros(size(sur.points,1),1);

for n = 1:size(sur.points,1)
    angleVectorToReferenceVector(n) = atan2(norm(cross(referenceVectorArray(n,:)',referenceVector')),dot(referenceVectorArray(n,:)',referenceVector'));
end

% find the vector which produces the smallest angle compared to the referenceVector
[value,idxPointOnSur] = min(angleVectorToReferenceVector);

referenceVectorMatch = referenceVectorArray(idxPointOnSur,:);

end
