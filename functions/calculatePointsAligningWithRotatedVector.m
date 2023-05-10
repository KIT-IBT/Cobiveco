function [pointRotated1, pointRotated2] = calculatePointsAligningWithRotatedVector(sur, referenceVector1, referenceVector2, projectedReferenceVectorArray)
% Calculate the vector closest aligning with referenceVector1 and reference_vector2
% defined best by minimizing the angle and define the point this vector points to.
%
% [pointRotated1, pointRotated2] = calculatePointsAligningWithRotatedVector(sur, referenceVector1, referenceVector2, projectedReferenceVectorArray))
%
% Inputs:
%
%   sur struct, surface structure representing the valve annuli with points being [numSurPoints x 3]
%   referenceVector1, double [1x3]
%   referenceVector2, double [1x3]
%   projectedReferenceVectorArray,
%
% Output:
%
%   Vector matching the referencevector best by minimizing the angle.
%   pointRotate1, double [1x3]
%   pointRotate2, double [1x3]
%
%
% Written by Lisa Pankewitz


% calculate set of vectors between center and each point in surface

angleVectorToReferenceVector1 = zeros(size(sur.points,1),1);

for n = 1:size(sur.points,1)
    angleVectorToReferenceVector1(n) = atan2(norm(cross(projectedReferenceVectorArray(n,:)',referenceVector1')),dot(projectedReferenceVectorArray(n,:)',referenceVector1'));
end

angleVectorToReferenceVector2 = zeros(size(sur.points,1),1);

for n = 1:size(sur.points,1)
    angleVectorToReferenceVector2(n) = atan2(norm(cross(projectedReferenceVectorArray(n,:)',referenceVector2')),dot(projectedReferenceVectorArray(n,:)',referenceVector2'));
end


[value,idxPointClosestToReferenceVector1] = min(abs(angleVectorToReferenceVector1));
pointRotated1 = sur.points(idxPointClosestToReferenceVector1,:);

[value,idxPointClosestToReferenceVector2] = min(abs(angleVectorToReferenceVector2));
pointRotated2 = sur.points(idxPointClosestToReferenceVector2,:);



end
