function plotValveAxis(sur, center, rotAxis, referenceVector, vectorRotated1, vectorRotated2, cfg, suffix)
% Plot a set of axis.
%
% plotValveAxis(sur, center, rotAxis, referenceVector, vectorRotated1, vectorRotated2, cfg, suffix)
%
% Inputs:
%   sur struct, surface representing the valve annuli with points being [numSurPoints x 3]
%   center, double [1x3]
%   rotAxis,double [1x3]
%   referenceVector, double [1x3]
%   vectorRotated1, double [1x3]
%   vectorRotated2, double [1x3]
%   cfg, struct, property of class cobiveco (see class description for details)
%   suffix, string
%
% Written by Lisa Pankewitz

vtkWrite(sur, sprintf('%ssur%s.vtp', cfg,suffix));
surAxes.points = [center; center+10*rotAxis; center+referenceVector; center+vectorRotated1; center+vectorRotated2];
surAxes.cells = int32([1 2; 1 3; 1 4; 1 5]);
surAxes.cellTypes = uint8([3; 3; 3; 3]);
surAxes.cellData.axis = uint8([1; 2; 3; 4]);
vtkWrite(surAxes, sprintf('%ssurAxes%s.vtp', cfg,suffix));

end