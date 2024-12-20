function M = cobiveco_computeMappingMatrix_briges(source, target, method, searchradius, verbose)
    % Uses ventricular coordinates to compute a matrix that maps from the 
    % points of a source mesh to the points of a target mesh.
    % Linear or nearest neighbor interpolation can be chosen for the mapping.
    %
    % Syntax:
    %  M = cobiveco_computeMappingMatrix(source, target, method, searchradius, verbose)
    %
    % Inputs:
    % - source:       VTK struct of source mesh containing cobiveco coordinates
    % - target:       VTK struct of target mesh containing cobiveco coordinates
    % - method:       interpolation method: 'linear' (default) or 'nearest'
    % - searchradius: radius used to search for source cell centroids;
    %                 unit: mean edge length; default: 2; 
    %                 only used for method=='linear'
    % - verbose:      whether to print status messages; default: false
    %
    % Output:
    % - M:            sparse mapping matrix (numTargetPoints x numSourcePoints)
    %
    % Written in 2020 by Steffen Schuler
    % Institute of Biomedical Engineering, KIT
    % www.ibt.kit.edu
    
    if nargin < 5
        verbose = false;
    end
    if nargin < 4 || isempty(searchradius)
        searchradius = 2; % unit: mean edge length
    end
    if nargin < 3 || isempty(method)
        method = 'linear';
    end
    
    source = updateOutliersWithNeighborValues(source);
    target = updateOutliersWithNeighborValues(target);

        
    %% Scale ventricular coords to have a similar change across one tet.
    
    if verbose, tic; fprintf('Scaling coordinates...                 '); end
    
    tv_cells = round(source.pointData.tv(source.cells));
    tv_norm = mean(vtkEdgeLengths(source)) / norm(max(source.points,[],1)-min(source.points,[],1));
    
    ab_cells = source.pointData.ab(source.cells);
    ab_norm = median(max(ab_cells,[],2) - min(ab_cells,[],2));

    rt_cells = source.pointData.rt(source.cells);
    rt_norm = median(max(rt_cells,[],2) - min(rt_cells,[],2));
    
    tm_cells = source.pointData.tm(source.cells);
    tm_norm = median(max(tm_cells,[],2) - min(tm_cells,[],2));
    tm_norm = max(tm_norm, 0.1);
    
    if verbose, fprintf('%.1f seconds\n', toc); end
    
    %%
    
    numSrcPoints = size(source.points,1);
    numTarPoints = size(target.pointData.ab,1);
    
    if strcmp(method, 'linear')
        %% Linear interpolation
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % For each target point, find the source cell centroid with the 
        % smallest euclidean distance in ventricular coords.
        
        if verbose, tic; fprintf('Searching closest centroids...         '); end
        
        X = NaN(size(source.cells,1), 4);
        X(:,1) = round(mean(tv_cells,2)) / tv_norm;
        X(:,2) = mean(ab_cells,2) / ab_norm;
        X(:,3) = mean(tm_cells,2) / tm_norm;
        % add for bridges, as for bridges no rotational will be used
        X(:,4) = mean(rt_cells,2) / rt_norm;
    
        Y = NaN(numTarPoints, 4);
        Y(:,1) = round(target.pointData.tv) / tv_norm;
        Y(:,2) = target.pointData.ab / ab_norm;
        Y(:,3) = target.pointData.tm / tm_norm;
        % adding rt for bridges where we do not use the cos and sin functions
        Y(:,4) = target.pointData.rt / rt_norm;
    
        Mdl1 = KDTreeSearcher(X);
        pointIds = knnsearch(Mdl1, Y);
    
        if verbose, fprintf('%.1f seconds\n', toc); end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % For each initial centroid, find all centroids within a radius
        % (euclidean distance in actual, cartesian coords).
        % This is done separately for left and right.

        % we might want to do this separately for bridges and not bridges, the criteria could be th apicobasal bigger than 1
    
        if verbose, tic; fprintf('Searching centroids within radius...   '); end
        
        s_cells = source.cells;
        s_points = double(source.points);
        s_centroids = squeeze(mean(reshape(s_points(s_cells,:),[],size(s_cells,2),size(s_points,2)),2));
    
        s_tv = round(mean(tv_cells,2));
        s_l = find(s_tv==0);
        s_r = find(s_tv==1);
    
        t_tv = round(target.pointData.tv);
        t_l = find(t_tv==0);
        t_r = find(t_tv==1);
    
        Mdl2l = KDTreeSearcher(s_centroids(s_l,:));
        Mdl2r = KDTreeSearcher(s_centroids(s_r,:));
    
        dist = searchradius * mean(vtkEdgeLengths(source));
        idx2l = rangesearch(Mdl2l, s_centroids(pointIds(t_l),:), dist);
        idx2r = rangesearch(Mdl2r, s_centroids(pointIds(t_r),:), dist);
    
        idx2 = cell(size(pointIds));
        for i = 1:numel(idx2l)
            idx2{t_l(i)} = s_l(idx2l{i});
        end
        for i = 1:numel(idx2r)
            idx2{t_r(i)} = s_r(idx2r{i});
        end
    
        if verbose, fprintf('%.1f seconds\n', toc); end
        if verbose, fprintf('%.1f seconds\n', toc); end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % For each target point, iterate over all source tets corresponding to 
        % the centroids and identify the tet to be used for interpolation.
        % For each candidate tet, barycentric coords reproducing the target
        % ventricular coords are computed. Bary coords are then used to 
        % identify the tet enclosing the point with the target coords or the
        % tet closest to this point.
        
        gcp; % start parallel pool, if not already running
        
        if verbose, tic; fprintf('Identifying cells for interpolation... '); end
    
        s_coords = [ ...
            source.pointData.ab ...
            source.pointData.tm ...
            source.pointData.rt ...
            ones(numSrcPoints,1) ...
            ]';
    
        t_coords = [ ...
            target.pointData.ab ...
            target.pointData.tm ...
            target.pointData.rt ...
            ones(numTarPoints,1) ...
            ]';
    
        numDims = size(source.cells,2);
        numCoords = size(s_coords,1);
        cellIds = NaN(numTarPoints,1);
        baryCoords = NaN(numTarPoints, numDims);
    
        baryMats = NaN(numDims, numCoords, size(s_cells,1));
        parfor i = 1:size(s_cells,1)
            A = s_coords(:,s_cells(i,:));
            baryMats(:,:,i) = pinv(A);
        end
        baryMats = permute(baryMats, [2 1 3]);
    
        parfor i = 1:numTarPoints
    
            candCellIds = idx2{i};
            if isempty(candCellIds)
                continue;
            end
    
            candBary = reshape(t_coords(:,i)' * reshape(baryMats(:,:,candCellIds), numCoords, []), numDims, [])';
    
            [~,k] = min(max(abs(candBary-0.5), [], 2));
    
            cellIds(i) = candCellIds(k);
            baryCoords(i,:) = candBary(k,:);
    
        end
    
        if verbose, fprintf('%.1f seconds\n', toc); end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Use barycentric coords to build the mapping matrix.
    
        if verbose, tic; fprintf('Building mapping matrix...             '); end
    
        nans = isnan(cellIds);
        cellIds(nans) = 1;
        baryCoords(nans,:) = NaN;
    
        i = reshape(repmat((1:numTarPoints)', 1, numDims), [], 1);
        j = double(reshape(s_cells(cellIds,:), [], 1));
        v = baryCoords(:);
        M = sparse(i, j, v, numTarPoints, numSrcPoints);
    
        if verbose, fprintf('%.1f seconds\n', toc); end
        
    elseif strcmp(method, 'nearest')
        %% Nearest neighbor interpolation
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % For each target point, find the source point with the smallest 
        % euclidean distance in ventricular coords.
        
        if verbose, tic; fprintf('Searching closest points...            '); end
        
        X = NaN(size(source.cells,1), 4);
        X(:,1) = round(mean(tv_cells,2)) / tv_norm;
        X(:,2) = mean(ab_cells,2) / ab_norm;
        X(:,3) = mean(tm_cells,2) / tm_norm;
        % add for bridges, as for bridges no rotational will be used
        X(:,4) = mean(rt_cells,2) / rt_norm;
    
        Y = NaN(numTarPoints, 4);
        Y(:,1) = round(target.pointData.tv) / tv_norm;
        Y(:,2) = target.pointData.ab / ab_norm;
        Y(:,3) = target.pointData.tm / tm_norm;
        % adding rt for bridges where we do not use the cos and sin functions
        Y(:,4) = target.pointData.rt / rt_norm;
    
        Mdl = KDTreeSearcher(X);
        pointIds = knnsearch(Mdl, Y);
    
        if verbose, fprintf('%.1f seconds\n', toc); end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Build mapping matrix.
    
        if verbose, tic; fprintf('Building mapping matrix...             '); end
    
        nans = isnan(pointIds);
        pointIds(nans) = 1;
    
        i = 1:numTarPoints;
        j = pointIds;
        v = ones(numTarPoints,1);
        v(nans) = NaN;
        M = sparse(i, j, v, numTarPoints, numSrcPoints);
    
        if verbose, fprintf('%.1f seconds\n', toc); end
        
    else
        error('Unknown method ''%s''', method);
    end
    
    %keyboard;
    
    end