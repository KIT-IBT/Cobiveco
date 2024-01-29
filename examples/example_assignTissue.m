%% example_assignTissue

% This script creates a tissue based on a CobivecoX coordinate
% in shape of an ellipsoid, and if given, several small spheres
% around the main ellipsoid.
% Note: The rotated result is needed to run this example.
%
% To define the tissue, the parameters in
% constructTissueEllipsoidWithSmallSpheres are needed.
%
% That is, the ellipse is defined by
%     centroidCobiveco = [tv, tm, ab, rt],
%     xSemiaxisLength, ySemiaxisLength, zSemiaxisLength
% and the spheres are defined by
%     numberOfSmallSpheres,
%     distanceBetweenCentroidAndSphereCenters
%     radiusSmallSphere
%     planeSpecificationForSpheres
%
% See also constructTissueEllipsoidWithSmallSpheres
%
% Simula 2022

% The  main idea is to
% label target points intersecting with an ellipsoid given as input as tissue.
% The ellipsoid may be surrounded by spheres which also should will labelled as tissue.

%% Define paths

addpath('result_example/');
sourcePrefix = 'result_example/';
outputPrefix = 'result_example/';

%% Create directory if it does not exist

outPath = fileparts(outputPrefix);
if ~isempty(outPath) && ~exist(outPath, 'dir')
    mkdir(outPath);
end

%% Load geometries

% Important that one uses the rotated result
% (heart axis parallel with xyz-axes in Cartesian)
source = vtkRead(sprintf('%sresultR.vtu', sourcePrefix));

%% Choose characteristics of the ellipsoid and small spheres

% Choosing centroid: centroidCobiveco = [tv, tm, ab, rt]
centroidCobiveco = [1, 0.3,0.845, 0.67];

% Choosing values for the semi-axes
[xSemiaxisLength, ySemiaxisLength, zSemiaxisLength] = deal(9,7,7);

% Choosing distanceBetweenCentroidAndSphereCenters
distanceBetweenCentroidAndSphereCenters = 5;

% Choosing numberOfSmallSpheres and radiusSmallSphere
numberOfSmallSpheres = 5; % aim for interger value (represented as double); 0 for non
radiusSmallSphere = 1.5;

planeSpecificationForSpheres = "yz";
%% Call constructTissueEllipsoid to flag tissue intersecting with the ellipsoid determined above

% Adds the flag "tissue" to all vertices in ellipsoid
tissueEllipsoidOnly = constructTissueEllipsoidWithoutSavingIndices(source,centroidCobiveco,xSemiaxisLength,ySemiaxisLength,zSemiaxisLength);

% Adds the flag "tissue" to all vertices in ellipsoid and small spheres
tissueIncludingSpheres = constructTissueEllipsoidWithSmallSpheres(source,centroidCobiveco,xSemiaxisLength,ySemiaxisLength,zSemiaxisLength,distanceBetweenCentroidAndSphereCenters,numberOfSmallSpheres,radiusSmallSphere,planeSpecificationForSpheres);

%% Export results
vtkWrite(tissueEllipsoidOnly, sprintf('%sresult_tissueEllipsoidOnly.vtu', outputPrefix));
vtkWrite(tissueIncludingSpheres, sprintf('%sresult_tissueIncludingSpheres.vtu', outputPrefix));

%% Helper function

% Probably not needed
function [xSemiaxisLength, ySemiaxisLength, zSemiaxisLength] = determingAppropriateLengthsOfSemiaxes(source)
% Determine appropriate lengths of semiaxes based on how big the given
% heart is.
%
% [xSemiaxisLength, ySemiaxisLength, zSemiaxisLength] = determingAppropriateLengthsOfSemiaxes(source)
%
% Inputs:
%   source: object provided with cobivecoX coordinates, structure
%
% Output:
%   xSemiaxisLength: length of semiaxes parallel to the x-axis, double
%   ySemiaxisLength: length of semiaxes parallel to the y-axis, double
%   zSemiaxisLength: length of semiaxes parallel to the z-axis, double
%
% Copyright 2022 Simula Research Laboratory

%     size(source.points) % nx3
    maxAbsX = max(abs(source.points(1,:)));
    maxAbsY = max(abs(source.points(2,:)));
    maxAbsZ = max(abs(source.points(3,:)));

    [xSemiaxisLength, ySemiaxisLength, zSemiaxisLength] = deal(0.15*maxAbsX,0.25*maxAbsY,0.15*maxAbsZ);
end
