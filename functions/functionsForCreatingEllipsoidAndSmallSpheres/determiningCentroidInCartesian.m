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
%
% Output:
%   centroidCartesian: closest point in mesh in Cartesian coordinates, double [1x3]   
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

    centroidCobiveco
    estimatedCentroidCobivecoCoords = double(estimatedCentroidCobivecoCoords) % convert from single to double

    centroidCartesian = source.points(index, :)
end