% This example script creates an AHA parcellation as described
% in M.D. Cerqueira et al., Circulation. 2002;105:539–542 (https://doi.org/10.1161/hc0402.102975)
% on a mesh that has cobiveco coordinates assigned to it.

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

%% Start AHA parcellation

for i=1:length(source.pointData.tv)
    if source.pointData.tv(i) == 0 %start with left ventricle
        if source.pointData.ab(i) > 7/10 % the basal region, divided in 6 sectors
            if source.pointData.rt(i) >= 1/2 && source.pointData.rt(i) < 2/3
                source.pointData.parcellation(i) = 1; %anterior
            elseif source.pointData.rt(i) >= 2/3 && source.pointData.rt(i) < 5/6
                source.pointData.parcellation(i) = 2; %anteroseptal
            elseif source.pointData.rt(i) >= 5/6 && source.pointData.rt(i) < 1
                source.pointData.parcellation(i) = 3; %inferoseptal
            elseif source.pointData.rt(i) == 1 || source.pointData.rt(i) < 1/6
                source.pointData.parcellation(i) = 4; %inferior
            elseif source.pointData.rt(i) >= 1/6 && source.pointData.rt(i) < 1/3
                source.pointData.parcellation(i) = 5; %inferolateral
            elseif source.pointData.rt(i) >= 1/3 && source.pointData.rt(i) < 1/2
                source.pointData.parcellation(i) = 6; %anterolateral
            end
        elseif source.pointData.ab(i) <= 7/10 && source.pointData.ab(i) > 4/10 % the mid region, divided in 6 sectors
            if source.pointData.rt(i) >= 1/2 && source.pointData.rt(i) < 2/3
                source.pointData.parcellation(i) = 7; %anterior
            elseif source.pointData.rt(i) >= 2/3 && source.pointData.rt(i) < 5/6
                source.pointData.parcellation(i) = 8; %anteroseptal
            elseif source.pointData.rt(i) >= 5/6 && source.pointData.rt(i) < 1
                source.pointData.parcellation(i) = 9; %inferoseptal
            elseif source.pointData.rt(i) == 1 || source.pointData.rt(i) < 1/6
                source.pointData.parcellation(i) = 10; %inferior
            elseif source.pointData.rt(i) >= 1/6 && source.pointData.rt(i) < 1/3
                source.pointData.parcellation(i) = 11; %inferolateral
            elseif source.pointData.rt(i) >= 1/3 && source.pointData.rt(i) < 1/2
                source.pointData.parcellation(i) = 12; %anterolateral
            end
        elseif source.pointData.ab(i) <= 4/10 && source.pointData.ab(i) > 1/10 % the apical region, divided in 4 sectors
             if source.pointData.rt(i) >= 11/24 && source.pointData.rt(i) < 17/24
                source.pointData.parcellation(i) = 13; % anterior
            elseif source.pointData.rt(i) >= 17/24 && source.pointData.rt(i) < 23/24
                source.pointData.parcellation(i) = 14; %septal
            elseif source.pointData.rt(i) >= 23/24 || source.pointData.rt(i) < 5/24
                source.pointData.parcellation(i) = 15; %inferior
            elseif source.pointData.rt(i) >= 5/24 && source.pointData.rt(i) < 11/24
                source.pointData.parcellation(i) = 16; %lateral
             end
        elseif source.pointData.ab(i) <= 1/10 %the apex region
            source.pointData.parcellation(i) = 17;
        end
    elseif source.pointData.tv(i) == 1 % right ventricle, same regions as in the left one
        if source.pointData.ab(i) > 7/10 % basal region, divided in 6 sectors
            if source.pointData.rt(i) >= 1/2 && source.pointData.rt(i) < 2/3
                source.pointData.parcellation(i) = 18; %anterior
            elseif source.pointData.rt(i) >= 2/3 && source.pointData.rt(i) < 5/6
                source.pointData.parcellation(i) = 19; %anteroseptal
            elseif source.pointData.rt(i) >= 5/6 && source.pointData.rt(i) < 1
                source.pointData.parcellation(i) = 20; %inferoseptal
            elseif source.pointData.rt(i) == 1 || source.pointData.rt(i) < 1/6
                source.pointData.parcellation(i) = 21; %inferior
            elseif source.pointData.rt(i) >= 1/6 && source.pointData.rt(i) < 1/3
                source.pointData.parcellation(i) = 22; %inferolateral
            elseif source.pointData.rt(i) >= 1/3 && source.pointData.rt(i) < 1/2
                source.pointData.parcellation(i) = 23; %anterolateral
            end
        elseif source.pointData.ab(i) <= 7/10 && source.pointData.ab(i) > 4/10 % the mid region, divided in 6 sectors
            if source.pointData.rt(i) >= 1/2 && source.pointData.rt(i) < 2/3
                source.pointData.parcellation(i) = 24; %anterior
            elseif source.pointData.rt(i) >= 2/3 && source.pointData.rt(i) < 5/6
                source.pointData.parcellation(i) = 25; %anteroseptal
            elseif source.pointData.rt(i) >= 5/6 && source.pointData.rt(i) < 1
                source.pointData.parcellation(i) = 26; %inferoseptal
            elseif source.pointData.rt(i) == 1 || source.pointData.rt(i) < 1/6
                source.pointData.parcellation(i) = 27; %inferior
            elseif source.pointData.rt(i) >= 1/6 && source.pointData.rt(i) < 1/3
                source.pointData.parcellation(i) = 28; %inferolateral
            elseif source.pointData.rt(i) >= 1/3 && source.pointData.rt(i) < 1/2
                source.pointData.parcellation(i) = 29; %anterolateral
            end
        elseif source.pointData.ab(i) <= 4/10 && source.pointData.ab(i) > 1/10 %the apical region, divided in 4 sectors
             if source.pointData.rt(i) >= 11/24 && source.pointData.rt(i) < 17/24
                source.pointData.parcellation(i) = 30; %anterior
            elseif source.pointData.rt(i) >= 17/24 && source.pointData.rt(i) < 23/24
                source.pointData.parcellation(i) = 31; %septal
            elseif source.pointData.rt(i) >= 23/24 || source.pointData.rt(i) < 5/24
                source.pointData.parcellation(i) = 32; %inferior
            elseif source.pointData.rt(i) >= 5/24 && source.pointData.rt(i) < 11/24
                source.pointData.parcellation(i) = 33; %lateral
             end
        elseif source.pointData.ab(i) <= 1/10 %the apex region
            source.pointData.parcellation(i) = 34;
        end
    end
end

%% Export result

source.pointData.parcellation = source.pointData.parcellation';
vtkWrite(source, sprintf('%sresultWithAHAPar_geo2.vtu', outputPrefix));
