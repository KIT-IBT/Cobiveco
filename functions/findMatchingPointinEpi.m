function [idxClosestPointsInEpicardium] = findMatchingPointinEpi(struct, arrayPathTvPv, nodesListEpicardium)

    % Finds the matching points defined in Endocardium in the epicardium
    %
    % [idxClosestPointsInEpicardium] = findMatchingPointinEpi(struct, arrayPath, nodesListEpicardium)
    % 
    % Inputs:
    %
    %   struct, [numSurPoints x 3]
    %   arrayPathTvPv, array double [nx1]
    %   nodesListEpicardium, double [nx1]
    %
    % Output:
    %   
    %    idxClosestPointsInEpicardium, double [nx1]
    %
    % Written by Lisa Pankewitz


    subsetPointsArrayPath1TvPv = struct.vol.points(struct.surToVol(arrayPathTvPv),:);
    [index, distance] = knnsearch(struct.vol.points(nodesListEpicardium,:),subsetPointsArrayPath1TvPv);
    idxClosestPointsInEpicardium = nodesListEpicardium(index);
end