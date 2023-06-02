% Example script creating a parcellation on a mesh
% with existing result(cobiveco assignment).
% The object is divided in for sectors: left anterior, left posterior, right posterior and right anterior.

addpath('result_geo2/');

%% Define paths

sourcePrefix = 'result_geo2/';
outputPrefix = 'result_parcellation/';

%% Create directory if it does not exist

outPath = fileparts(outputPrefix);
if ~isempty(outPath) && ~exist(outPath, 'dir')
    mkdir(outPath);
end

%% Load geometries

source = vtkRead(sprintf('%sresult.vtu', sourcePrefix));

%% Start parcellation

for i=1:length(source.pointData.tv)
    if source.pointData.tv(i) == 0 && (0 <= source.pointData.rt(i) && source.pointData.rt(i) <= 1/3 || 5/6 < source.pointData.rt(i) && source.pointData.rt(i) <= 1)
        source.pointData.parcellation(i) = 1; %left anterior
    elseif source.pointData.tv(i) == 0 &&  1/3 < source.pointData.rt(i) && source.pointData.rt(i) <= 5/6
        source.pointData.parcellation(i) = 2; %left posterior
    elseif source.pointData.tv(i) == 1 &&  1/3 < source.pointData.rt(i) && source.pointData.rt(i) <= 5/6
        source.pointData.parcellation(i) = 3; %right posterior
    elseif source.pointData.tv(i) == 1 && (0 <= source.pointData.rt(i) && source.pointData.rt(i) <= 1/3 || 5/6 < source.pointData.rt(i) && source.pointData.rt(i) <= 1)
        source.pointData.parcellation(i) = 4; %right anterior
    end
end

%% Export result

source.pointData.parcellation = source.pointData.parcellation';
vtkWrite(source, sprintf('%sresultWithPar_geo2.vtu', outputPrefix));