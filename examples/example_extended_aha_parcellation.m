% This example script serves the purpose of creating an AHA parcellation as described
% in M.D. Cerqueira et al., Circulation. 2002;105:539–542 (https://doi.org/10.1161/hc0402.102975)
% on a result previously found thanks to cobivecoX.
% In script two new regions has been added to the ones from the
% previous model: the bridges and the so called "upperbasal" region.

addpath('result_templateHLHS/');

%% Define paths

sourcePrefix = 'result_templateHLHS/';
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
        if source.pointData.ab(i) > 1 %the bridge region, divided in 2
            if source.pointData.ab(i) >= 1.2
                source.pointData.parcellation(i) = 1; % anterior
            else
                source.pointData.parcellation(i) = 2; %inferior
            end
        elseif source.pointData.ab(i) > 4/5 && source.pointData.ab(i) <= 1 %the "upperbasal" region, divided in 6 sectors
            if source.pointData.rt(i) >= 87/200 && source.pointData.rt(i) < 29/50
                source.pointData.parcellation(i) = 3; %anterior
            elseif source.pointData.rt(i) >= 29/50 && source.pointData.rt(i) < 71/100
                source.pointData.parcellation(i) = 4; %anteroseptal
            elseif source.pointData.rt(i) >= 71/100 && source.pointData.rt(i) < 1
                source.pointData.parcellation(i) = 5; %inferoseptal
            elseif source.pointData.rt(i) == 1 || source.pointData.rt(i) < 29/200
                source.pointData.parcellation(i) = 6; %inferior
            elseif source.pointData.rt(i) >= 29/200 && source.pointData.rt(i) < 29/100
                source.pointData.parcellation(i) = 7; %inferolateral
            elseif source.pointData.rt(i) >= 29/100 && source.pointData.rt(i) < 87/200
                source.pointData.parcellation(i) = 8; %anterolateral
            end
        elseif source.pointData.ab(i) > 14/25 && source.pointData.ab(i) <= 4/5 %the basal region, divided in 6 sectors
            if source.pointData.rt(i) >= 87/200 && source.pointData.rt(i) < 29/50
                source.pointData.parcellation(i) = 9; %anterior
            elseif source.pointData.rt(i) >= 29/50 && source.pointData.rt(i) < 71/100
                source.pointData.parcellation(i) = 10; %anteroseptal
            elseif source.pointData.rt(i) >= 71/100 && source.pointData.rt(i) < 1
                source.pointData.parcellation(i) = 11; %inferoseptal
            elseif source.pointData.rt(i) == 1 || source.pointData.rt(i) < 29/200
                source.pointData.parcellation(i) = 12; %inferior
            elseif source.pointData.rt(i) >= 29/200 && source.pointData.rt(i) < 29/100
                source.pointData.parcellation(i) = 13; %inferolateral
            elseif source.pointData.rt(i) >= 29/100 && source.pointData.rt(i) < 87/200
                source.pointData.parcellation(i) = 14; %anterolateral
            end
        elseif source.pointData.ab(i) <= 14/25 && source.pointData.ab(i) > 8/25 %the mid region, divided in 6 sectors
            if source.pointData.rt(i) >= 87/200 && source.pointData.rt(i) < 29/50
                source.pointData.parcellation(i) = 15; %anterior
            elseif source.pointData.rt(i) >= 29/50 && source.pointData.rt(i) < 71/100
                source.pointData.parcellation(i) = 16; %anteroseptal
            elseif source.pointData.rt(i) >= 71/100 && source.pointData.rt(i) < 1
                source.pointData.parcellation(i) = 17; %inferoseptal
            elseif source.pointData.rt(i) == 1 || source.pointData.rt(i) < 29/200
                source.pointData.parcellation(i) = 18; %inferior
            elseif source.pointData.rt(i) >= 29/200 && source.pointData.rt(i) < 29/100
                source.pointData.parcellation(i) = 19; %inferolateral
            elseif source.pointData.rt(i) >= 29/100 && source.pointData.rt(i) < 87/200
                source.pointData.parcellation(i) = 20; %anterolateral
            end
        elseif source.pointData.ab(i) <= 8/25 && source.pointData.ab(i) > 2/25 %the apical region, divided in 4 sectors
            if source.pointData.rt(i) >= 83/200 && source.pointData.rt(i) < 133/200
                source.pointData.parcellation(i) = 21; %anterior
            elseif source.pointData.rt(i) >= 133/200 && source.pointData.rt(i) < 183/200
                source.pointData.parcellation(i) = 22; %septal
            elseif source.pointData.rt(i) >= 183/200 || source.pointData.rt(i) < 33/200
                source.pointData.parcellation(i) = 23; %inferior
            elseif source.pointData.rt(i) >= 33/200 && source.pointData.rt(i) < 83/200
                source.pointData.parcellation(i) = 24; %lateral
             end
        elseif source.pointData.ab(i) <= 2/25 %the apex region
            source.pointData.parcellation(i) = 25;
        end
    elseif source.pointData.tv(i) == 1 %right ventricle, same regions as in the left one
        if source.pointData.ab(i) > 1 %the bridge region, divided in 2
            if 0.45 < source.pointData.rt(i) && source.pointData.rt(i) <= 1
                source.pointData.parcellation(i) = 26; %anterior
            elseif -0.1 <= source.pointData.rt(i) && source.pointData.rt(i) <= 0.45
                source.pointData.parcellation(i) = 27; %inferior
            end
        elseif source.pointData.ab(i) > 4/5 && source.pointData.ab(i) <= 1 %the "upperbasal" region, divided in 6 sectors
            if source.pointData.rt(i) >= 87/200 && source.pointData.rt(i) < 29/50
                source.pointData.parcellation(i) = 28; %anterior
            elseif source.pointData.rt(i) >= 29/50 && source.pointData.rt(i) < 71/100
                source.pointData.parcellation(i) = 29; %anteroseptal
            elseif source.pointData.rt(i) >= 71/100 && source.pointData.rt(i) < 1
                source.pointData.parcellation(i) = 30; %inferoseptal
            elseif source.pointData.rt(i) == 1 || source.pointData.rt(i) < 29/200
                source.pointData.parcellation(i) = 31; %inferior
            elseif source.pointData.rt(i) >= 29/200 && source.pointData.rt(i) < 29/100
                source.pointData.parcellation(i) = 32; %inferolateral
            elseif source.pointData.rt(i) >= 29/100 && source.pointData.rt(i) < 87/200
                source.pointData.parcellation(i) = 33; %anterolateral
            end
        elseif source.pointData.ab(i) > 14/25 && source.pointData.ab(i) <= 4/5 %the basal region, divided in 6 sectors
            if source.pointData.rt(i) >= 87/200 && source.pointData.rt(i) < 29/50
                source.pointData.parcellation(i) = 34; %anterior
            elseif source.pointData.rt(i) >= 29/50 && source.pointData.rt(i) < 71/100
                source.pointData.parcellation(i) = 35; %anteroseptal
            elseif source.pointData.rt(i) >= 71/100 && source.pointData.rt(i) < 1
                source.pointData.parcellation(i) = 36; %inferoseptal
            elseif source.pointData.rt(i) == 1 || source.pointData.rt(i) < 29/200
                source.pointData.parcellation(i) = 37; %inferior
            elseif source.pointData.rt(i) >= 29/200 && source.pointData.rt(i) < 29/100
                source.pointData.parcellation(i) = 38; %inferolateral
            elseif source.pointData.rt(i) >= 29/100 && source.pointData.rt(i) < 87/200
                source.pointData.parcellation(i) = 39; %anterolateral
            end
        elseif source.pointData.ab(i) <= 14/25 && source.pointData.ab(i) > 8/25 %the mid region, divided in 6 sectors
            if source.pointData.rt(i) >= 87/200 && source.pointData.rt(i) < 29/50
                source.pointData.parcellation(i) = 40; %anterior
            elseif source.pointData.rt(i) >= 29/50 && source.pointData.rt(i) < 71/100
                source.pointData.parcellation(i) = 41; %anteroseptal
            elseif source.pointData.rt(i) >= 71/100 && source.pointData.rt(i) < 1
                source.pointData.parcellation(i) = 42; %inferoseptal
            elseif source.pointData.rt(i) == 1 || source.pointData.rt(i) < 29/200
                source.pointData.parcellation(i) = 43; %inferior
            elseif source.pointData.rt(i) >= 29/200 && source.pointData.rt(i) < 29/100
                source.pointData.parcellation(i) = 44; %inferolateral
            elseif source.pointData.rt(i) >= 29/100 && source.pointData.rt(i) < 87/200
                source.pointData.parcellation(i) = 45; %anterolateral
            end
        elseif source.pointData.ab(i) <= 8/25 && source.pointData.ab(i) > 2/25 %the apical region, divided in 4 sectors
            if source.pointData.rt(i) >= 83/200 && source.pointData.rt(i) < 133/200
                source.pointData.parcellation(i) = 46; %anterior
            elseif source.pointData.rt(i) >= 133/200 && source.pointData.rt(i) < 183/200
                source.pointData.parcellation(i) = 47; %septal
            elseif source.pointData.rt(i) >= 183/200 || source.pointData.rt(i) < 33/200
                source.pointData.parcellation(i) = 48; %inferior
            elseif source.pointData.rt(i) >= 33/200 && source.pointData.rt(i) < 83/200
                source.pointData.parcellation(i) = 49; %lateral
             end
        elseif source.pointData.ab(i) <= 2/25 %the apex region
            source.pointData.parcellation(i) = 50;
        end
    end
end

%% Export result

source.pointData.parcellation = source.pointData.parcellation';
vtkWrite(source, sprintf('%sresultWithAHAPar_templateHLHS.vtu', outputPrefix));
