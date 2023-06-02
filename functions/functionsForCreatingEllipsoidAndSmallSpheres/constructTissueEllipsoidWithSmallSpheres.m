function source = constructTissueEllipsoidWithSmallSpheres(source,centroidCobiveco,xSemiaxisLength,ySemiaxisLength,zSemiaxisLength,distanceBetweenCentroidAndSphereCenters,numberOfSmallSpheres,radiusSmallSphere,planeSpecificationForSpheres)
% Adds the flag "tissue" to all vertices close to the centroid in the form
% of an ellipsoid with axes lengths 
%       xSemiaxisLength, ySemiaxisLength and zSemiaxisLength
% and spheres described by input.
% That is, a given number of spheres with centeres (in the XY-plane) with given distance to 
% the centroid and given radius.
%
% source = constructTissueEllipsoidWithSmallSpheres(source,centroidCobiveco,xSemiaxisLength,ySemiaxisLength,zSemiaxisLength,distanceBetweenCentroidAndSphereCenters,numberOfSmallSpheres,radiusSmallSphere)
%
% Inputs:
%   source: object provided with cobivecoX coordinates, structure
%   centroidCobiveco: center of the ellipsoid in Cobiveco coordinates, double [1x4]
%   xSemiaxisLength: chosen elipsoid x-semiaxis length, double
%   ySemiaxisLength: chosen elipsoid y-semiaxis length, double
%   zSemiaxisLength: chosen elipsoid z-semiaxis length, double
%   distanceBetweenCentroidAndSphereCenters
%   numberOfSmallSpheres:
%   radiusSmallSphere: double
%
% Output:
%   source: object with new flag "tissue"
% 
% See also determiningCentroidInCartesian, 
%   determiningSphereCentersCartesian, 
%   determiningPointIndicesTissueEllipsoid
%
% Simula 2022

% Need to ensure that numberOfSmallSpheres is non-negative and prefereably
% an int (but the latter lead to error if numberOfSmallSpheres = int(4))

%% Find centroid in Cartesian coords 
% since equation for ellipsoids are known in Cartesian but not in cobiveco.
centroidCartesian = determiningCentroidInCartesian(source,centroidCobiveco);
if isempty(centroidCartesian) % Just to ensure centroid is in Catesian coords, choose the origo        
    centroidCartesian = [0,0,0];
end

%% Ellipsoid
% Initialize pointIndicesForTissue by determining (point indices for) ellipsoid
pointIndicesForTissue = determiningPointIndicesTissueEllipsoid(source,centroidCartesian,xSemiaxisLength,ySemiaxisLength,zSemiaxisLength);
% size(pointIndicesForTissue)

%% Small spheres
if (numberOfSmallSpheres > 0) % Actually negative number would work, but need to ensure the number is not rounded to 0
    sphereCentersCartesian = determiningSphereCentersCartesian(source,centroidCartesian,distanceBetweenCentroidAndSphereCenters,numberOfSmallSpheres,planeSpecificationForSpheres)

    % Update pointIndicesForTissue by determining (point indices for) the small spheres
    for i = 1:numberOfSmallSpheres
         pointIndicesSphere = determiningPointIndicesTissueEllipsoid(source,sphereCentersCartesian(i,:),radiusSmallSphere,radiusSmallSphere,radiusSmallSphere);
         pointIndicesForTissue = [pointIndicesForTissue, pointIndicesSphere];
    %      size(pointIndicesForTissue)
    end

end

%% Mark tissue (Inefficient way of marking tissue)
% First mark tissue as 0 at all points, 
% then update indices found when calling
% determiningPointIndicesTissueEllipsoid, 
% ie., all points which are part of ellipsoid and small spheres

numberOfDataPoints = size(source.points, 1);
for i = 1:numberOfDataPoints
    source.pointData.tissue(i,1) = 0;
end 

for i = 1:size(pointIndicesForTissue, 2) 
    source.pointData.tissue(pointIndicesForTissue(1,i),1) = 1;
end

end
