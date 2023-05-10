addpath('result_geo2/');

sourcePrefix = 'result_geo2/';
outputPrefix = 'result_parcellation/';

outPath = fileparts(outputPrefix);
if ~isempty(outPath) && ~exist(outPath, 'dir')
    mkdir(outPath);
end

source = vtkRead(sprintf('%sresult.vtu', sourcePrefix));

for i=1:length(source.pointData.tv)
    if source.pointData.tv(i) == 0 && (0 <= source.pointData.rt(i) && source.pointData.rt(i) <= 1/3 || 5/6 < source.pointData.rt(i) && source.pointData.rt(i) <= 1)
        source.pointData.parcellation(i) = 1;
    elseif source.pointData.tv(i) == 0 &&  1/3 < source.pointData.rt(i) && source.pointData.rt(i) <= 5/6
        source.pointData.parcellation(i) = 2;
    elseif source.pointData.tv(i) == 1 &&  1/3 < source.pointData.rt(i) && source.pointData.rt(i) <= 5/6
        source.pointData.parcellation(i) = 3;
    elseif source.pointData.tv(i) == 1 && (0 <= source.pointData.rt(i) && source.pointData.rt(i) <= 1/3 || 5/6 < source.pointData.rt(i) && source.pointData.rt(i) <= 1)
        source.pointData.parcellation(i) = 4;
    end
end
source.pointData.parcellation = source.pointData.parcellation';
%par= source.pointData.parcellation;
vtkWrite(source, sprintf('%sresultWithPar_geo2.vtu', outputPrefix));