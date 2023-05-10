function source = createTissueArea(source,centroid, R, a, b, c)
% Adds the field "tissue" to all vertices close to the centroid in the form
% of an ellipsoid with axes length a, b and c.
%
% source = createTissueArea(source,centroid, R, a, b, c)
%
% Inputs:
%   source(struct): object provided with cobivecoX coordinates
%   centroid(double, [1x4]): center of the elipsoid
%   R(double, [4x4]): rotation matrix to fit the heart axes
%   a (int): chosen elipsoid axis length 
%   b (int): chosen elipsoid axis length
%   c (int): chosen elipsoid axis length
%
% Output:
%   source(struct): object with new field "tissue"
% written 2022, Simula Research Laboratory

distances_a = zeros(length(source.points),1);
distances_b = zeros(length(source.points),1);
distances_c = zeros(length(source.points),1);

% Rotate to fit the heart axes
coordinates = R * [source.points, zeros(length(source.points),1)]';
centroid = R * centroid';

% Equation for ellipsoid: ((x-x_c)/a)^2+((y-y_c)/b)^2+((z-z_c)/c)^2 <= 1
for i=1:length(coordinates')
    distances_a(i) = ((coordinates(1,i)-centroid(1))/a)^2;
    distances_b(i) = ((coordinates(2,i)-centroid(2))/b)^2;
    distances_c(i) = ((coordinates(3,i)-centroid(3))/c)^2;
    if distances_a(i)+distances_b(i)+distances_c(i) <= 1
        source.pointData.tissue(i,1) = 1;
    else 
        source.pointData.tissue(i,1) = 0;
    end
end
end