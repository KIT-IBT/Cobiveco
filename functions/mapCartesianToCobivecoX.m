function outputCobivecoX = mapCartesianToCobivecoX(inputStruct)
% Maps the cartesian coordinates of the given struct to the CobivecoX coordinates.
%
% outputCobivecoX = mapCartesianToCobivecoX(inputStruct)
%
% Input:
%   inputStruct: instance of the class cobiveco.
%
% Output:
%   outputCobivecoX (double): matrix with cobivecoX coordinates.

outputCobivecoX = zeros(length(inputStruct.points),4);
for i=1:length(inputStruct.points)
    outputCobivecoX(i,1) = inputStruct.pointData.tv(i); 
    outputCobivecoX(i,2) = inputStruct.pointData.tm(i);
    outputCobivecoX(i,3) = inputStruct.pointData.rt(i);
    outputCobivecoX(i,4) = inputStruct.pointData.ab(i); 
end
end
