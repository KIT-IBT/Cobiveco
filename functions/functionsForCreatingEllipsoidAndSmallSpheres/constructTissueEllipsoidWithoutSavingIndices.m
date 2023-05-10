function source = constructTissueEllipsoidWithoutSavingIndices(source,centroidCobiveco,xSemiaxisLength,ySemiaxisLength,zSemiaxisLength)

% Adds the flag "tissue" to all vertices close to the centroid in the form
% of an ellipsoid with axes lengths 
%       xSemiaxisLength, ySemiaxisLength and zSemiaxisLength
%
% source = constructTissueEllipsoidWithoutSavingIndices(source,centroidCobiveco,xSemiaxisLength,ySemiaxisLength,zSemiaxisLength)
%
% Inputs:
%   source: object provided with cobivecoX coordinates, structure
%   centroidCobiveco: center of the ellipsoid in Cobiveco coordinates, double [1x4]
%   xSemiaxisLength: chosen elipsoid x-semiaxis length, double
%   ySemiaxisLength: chosen elipsoid y-semiaxis length, double
%   zSemiaxisLength: chosen elipsoid z-semiaxis length, double
%
% Output:
%   source: object with new flag "tissue"
%
% Note that when xSemiaxisLength = ySemiaxisLength = zSemiaxisLength,
% results in a sphere with radius xSemiaxisLength.
%
% Simula 2022

% Let source.pointData.tissue be an nx1-array. 
% source.points is nx3 (Cartesian)
    numberOfDataPoints = size(source.points, 1);

% Find centroid in Cartesian coords 
% since equation for ellipsoids are known in Cartesian but not in
% cobiveco.
centroidCartesian = determiningCentroidInCartesian(source,centroidCobiveco);
if isempty(centroidCartesian) % Just to ensure centroid is in Catesian coords, choose the origo        
    centroidCartesian = [0,0,0];
end

% Define ellipsoid:
% Let  a = x_semiaxis_length, b = y_semiaxis_length, c = z_semiaxis_length. 
% Equation for ellipsoid: ((x-x_c)/a)^2+((y-y_c)/b)^2+((z-z_c)/c)^2 <= 1
for i = 1:numberOfDataPoints
    distanceXdir = ((source.points(i,1) - centroidCartesian(1))/xSemiaxisLength)^2;
    distanceYdir = ((source.points(i,2) - centroidCartesian(2))/ySemiaxisLength)^2;
    distanceZdir = ((source.points(i,3) - centroidCartesian(3))/zSemiaxisLength)^2;

    if distanceXdir + distanceYdir + distanceZdir <= 1
        source.pointData.tissue(i,1) = 1;
    else 
        source.pointData.tissue(i,1) = 0;
    end
end

end


%% Helper functions

function centroidCartesian = determiningCentroidInCartesian(source,centroidCobiveco)
% Given centroid in Cobiveco coordinates, this function minimizes the
% distance between given centroid input and a point existising in (the given)
% covibeco object. The closest point is returned in Cartesian coordinates.
%
% Here "closest" is the vertex in the meshmesh where the difference between
% wanted/input coordinates and existing coordinate is minimized by
%   min(sum(differences)), where differences stores the absolute values of
%   difference between two points.
%
% centroidCartesian = determiningCentroidInCartesian(source,centroidCobiveco)
%
% Inputs:
%   source: object provided with cobivecoX coordinates, struct
%   centroid: center of the elipsoid in Cobiveco coordinates, double [1x4]

% Output:
%   centroidCartesian: closest point in mesh in Cartesian coordinates   
%
% Simula 2022

%     % source.points is nx3 (Cartesian)
%     numberOfDataPoints = size(source.points, 1);

    % Assume index i in source.points 
    % is mapped to the same index in source.pointData
    % Locate index/indices to point closest to input centroid
    % Extract corresponding index/indices in source.points
     
    % Calculating absolute values of distances betweeen points in CobivecoX
    % coords to minimize them. This is wanted because the index is needed 
    % to deterine Cartesian coords of the centroid
    differences(:, 1) = source.pointData.tv - centroidCobiveco(1,1);
    differences(:, 2) = source.pointData.tm - centroidCobiveco(1,2);
    differences(:, 3) = source.pointData.ab - centroidCobiveco(1,3);
    differences(:, 4) = source.pointData.rt - centroidCobiveco(1,4);
    size(differences); % nx4
    differences = abs(differences);
    [minSumDifferences, index] = min(sum(differences,2)); % minSumDifferences is single not double
   
    estimatedCentroidCobivecoCoords(1) = source.pointData.tv(index);
    estimatedCentroidCobivecoCoords(2) = source.pointData.tm(index);
    estimatedCentroidCobivecoCoords(3) = source.pointData.ab(index);
    estimatedCentroidCobivecoCoords(4) = source.pointData.rt(index);

%     centroidCobiveco
    estimatedCentroidCobivecoCoords = double(estimatedCentroidCobivecoCoords); % convert from single to double

    centroidCartesian = source.points(index, :);
end