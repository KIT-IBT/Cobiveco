function sphereCentersCartesian = determiningSphereCentersCartesian(source,centroidCartesian,distanceBetweenCentroidAndSphereCenters,numberOfSmallSpheres,planeSpecificationForSpheres)
% For a given centroid in Cartesian coords, 
% this function determines the centers in Cartesian coords 
% for a given number of small spheres.
% These centers are either in the xy-, xz- or yz-plane (the former by
% default).
%
% source,centroidCartesian,distanceBetweenCentroidAndSphereCenters,numberOfSmallSpheres,planeSpecification)
%
% Inputs:
%   source: object provided with cobivecoX coordinates, structure
%   centroidCatesian: center of the ellipsoid in Cartesian coordinates, double [1x3]
%   numberOfSmallSpheres:
%   radiusSmallSphere: double
%   planeSpecificationForSpheres: one of the following strings specifying 
%           plane the centers are in: "xy", "xz", "yz" (xy is default)
%
% Output:
%   sphereCentersCartesian: Cartesian coordinates of spheres, [numberOfSmallSpheres x 3]
%  
% Simula 2022

    % This function is not called unless numberOfSmallSpheres > 0. 
    % It should also work for numberOfSmallSpheres ~= 0
    
    % Using polar coords in the xy-plane
    r = distanceBetweenCentroidAndSphereCenters;
    % Use angles in radians
    phi = 2*pi/numberOfSmallSpheres;

    if planeSpecificationForSpheres == "xz" 
        for i = 1:numberOfSmallSpheres
            sphereCentersCartesian(i,:) = ...
                [centroidCartesian(1) + r*cos(phi*(i-1)), ...
                 centroidCartesian(2), ...
                 centroidCartesian(3) + r*sin(phi*(i-1))];
        end

    elseif planeSpecificationForSpheres == "yz" 
        for i = 1:numberOfSmallSpheres
            sphereCentersCartesian(i,:) = ...
                [centroidCartesian(1), ...
                 centroidCartesian(2) + r*cos(phi*(i-1)), ...
                 centroidCartesian(3) + r*sin(phi*(i-1))];
        end

    else % xy-plane
        for i = 1:numberOfSmallSpheres
            sphereCentersCartesian(i,:) = ...
                [centroidCartesian(1) + r*cos(phi*(i-1)), ...
                 centroidCartesian(2) + r*sin(phi*(i-1)), ...
                 centroidCartesian(3)];
        end
    end

end