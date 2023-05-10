% This example script serves the purpose of creating a parcellation on a
% result previously found thanks to cobiveco. 
% The object is divided in for sectors: left anterior, left posterior, right posterior and right anterior.

addpath('result_CHD0017001/');
%% Define paths

sourcePrefix = 'result_CHD0017001/';
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
    if source.pointData.ab(i) <= 1 %start with ventricles
        if source.pointData.tv(i) == 0 && (0 <= source.pointData.rt(i) && source.pointData.rt(i) <= 1/3 || 2/3 < source.pointData.rt(i) && source.pointData.rt(i) <= 1)
            source.pointData.parcellation(i) = 1; %left anterior
        elseif source.pointData.tv(i) == 0 &&  1/3 < source.pointData.rt(i) && source.pointData.rt(i) <= 2/3
            source.pointData.parcellation(i) = 2; %left posterior
        elseif source.pointData.tv(i) == 1 &&  1/3 < source.pointData.rt(i) && source.pointData.rt(i) <= 2/3
            source.pointData.parcellation(i) = 3; %right anterior
        elseif source.pointData.tv(i) == 1 && (0 <= source.pointData.rt(i) && source.pointData.rt(i) <= 1/3 || 2/3 < source.pointData.rt(i) && source.pointData.rt(i) <= 1)
            source.pointData.parcellation(i) = 4; %right posterior
        end
    else %in the bridges
        if source.pointData.tv(i) == 1 && 0.45 < source.pointData.rt(i) && source.pointData.rt(i) <= 1
            source.pointData.parcellation(i) = 3; %right anterior
        elseif source.pointData.tv(i) == 1 && -0.1 <= source.pointData.rt(i) && source.pointData.rt(i) <= 0.45
            source.pointData.parcellation(i) = 4; %right posterior
        elseif source.pointData.ab(i) >= 1.2 
            source.pointData.parcellation(i) = 1; %left anterior
        else 
            source.pointData.parcellation(i) = 2; %left posterior
        end
    end
end 

%% Export result

source.pointData.parcellation = source.pointData.parcellation';
vtkWrite(source, sprintf('%sresultWithPar.vtu', outputPrefix));
