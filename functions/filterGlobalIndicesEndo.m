function [idxRv, idxLv, idxOppositeRv, idxOppositeLv] = filterGlobalIndicesEndo(struct1, idxArray, idsSept)

    % Define indices of LV and RV Nodes excluding the Septum to be excluded in
    % the (Sub)Graph Search
    %
    % [idxRv, idxLv, idxOppositeRv, idxOppositeLv] =
    % filterGlobalIndicesEndo(struct1, idxArray, idsSept)
    %
    % Inputs:
    %   struct1 [numPoints x 3]
    %   idxArray global volume indices of nodes in graph, int [numPoints x 1]
    %   idsSept arrays of septal IDs, int[numPoints x 1]
    %
    % Output:
    %   idxRv
    %   idxLv
    %   idxOppositeRv
    %   idxOppositeLv
    %
    % Written by Lisa Pankewitz
    
    % Node Filter IDs Epicardium
    nodesEpi = struct1.surToVol(struct1.sur.pointData.class == 1 | struct1.sur.pointData.class == 2);

    % Define Node IDs RV
    idsRvVol = find(struct1.tv == 1);

    % Exclude Septum from ID selection
    idsRvVol = setdiff(idsRvVol,idsSept);

    % Exclude Epicardium
    idsRvVol = setdiff(idsRvVol,nodesEpi);

    % Find Nodes common in Map and RV subsection
    [commonValuesRv, idxRv] = intersect(idxArray,idsRvVol,'stable');
    [commonValuesLv, idxOppositeRv] = setxor(idxArray,idsRvVol,'stable');

    % Repeat Process for LV Nodes
    % Define Node IDs LV
    idsLvVol = find(struct1.tv == 0);
    idsLvVol = setdiff(idsLvVol,idsSept);
    idsLvVol = setdiff(idsLvVol,nodesEpi);
    [commonValuesLv, idxLv] = intersect(idxArray,idsLvVol,'stable');
    [commonValuesLv, idxOppositeLv] = setxor(idxArray,idsLvVol,'stable');
 
    end
    

