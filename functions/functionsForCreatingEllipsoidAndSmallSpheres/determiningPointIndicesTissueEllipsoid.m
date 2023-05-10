function pointIndicesTissueEllipsoid = determiningPointIndicesTissueEllipsoid(source,centroidCartesian,xSemiaxisLength,ySemiaxisLength,zSemiaxisLength)
% Collects indices of all vertices close to the centroid in the form
% of an ellipsoid with axes lengths 
%       xSemiaxisLength, ySemiaxisLength and zSemiaxisLength
%
% pointIndicesTissueEllipsoid = constructTissueEllipsoid(source,centroidCobiveco,xSemiaxisLength,ySemiaxisLength,zSemiaxisLength)
%
% Inputs:
%   source: object provided with cobivecoX coordinates, structure
%   centroidCobiveco: center of the ellipsoid in Cobiveco coordinates, double [1x4]
%   xSemiaxisLength: chosen elipsoid x-semiaxis length, double
%   ySemiaxisLength: chosen elipsoid y-semiaxis length, double
%   zSemiaxisLength: chosen elipsoid z-semiaxis length, double
%
% Output:
%   pointIndicesTissueEllipsoid: indices which are to be flaged later,
%   array of integers
%
% Note that when xSemiaxisLength = ySemiaxisLength = zSemiaxisLength,
% results in a sphere with radius xSemiaxisLength.
%
% Simula 2022

% Let source.pointData.tissue be an nx1-array. 
% source.points is nx3 (Cartesian)
numberOfDataPoints = size(source.points, 1);

% Define ellipsoid:
% Let  a = x_semiaxis_length, b = y_semiaxis_length, c = z_semiaxis_length. 
% Equation for ellipsoid: ((x-x_c)/a)^2+((y-y_c)/b)^2+((z-z_c)/c)^2 <= 1
pointIndicesTissueEllipsoid = []; % Initiating array of indices we want to keep
for i = 1:numberOfDataPoints
    distanceXdir = ((source.points(i,1) - centroidCartesian(1))/xSemiaxisLength)^2;
    distanceYdir = ((source.points(i,2) - centroidCartesian(2))/ySemiaxisLength)^2;
    distanceZdir = ((source.points(i,3) - centroidCartesian(3))/zSemiaxisLength)^2;

    if distanceXdir + distanceYdir + distanceZdir <= 1
        pointIndicesTissueEllipsoid = [pointIndicesTissueEllipsoid, i];
    end
end
pointIndicesTissueEllipsoid = int32(pointIndicesTissueEllipsoid);

end