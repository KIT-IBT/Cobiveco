function M = cobiveco_computeMappingMatrix(source, target, method, searchradius, verbose)
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

%% Transform rt into rtSin and rtCos, if necessary

if ~isfield(source.pointData, 'rtSin') || ~isfield(source.pointData, 'rtCos')
    source.pointData.rt = min(max(source.pointData.rt,0),1);
    source.pointData.rtSin = sin(2*pi*source.pointData.rt);
    source.pointData.rtCos = cos(2*pi*source.pointData.rt);
end
if ~isfield(target.pointData, 'rtSin') || ~isfield(target.pointData, 'rtCos')
    target.pointData.rt = min(max(target.pointData.rt,0),1);
    target.pointData.rtSin = sin(2*pi*target.pointData.rt);
    target.pointData.rtCos = cos(2*pi*target.pointData.rt);
end

% Process source.pointData.ab
sourceAb = source.pointData.ab;
sourceRt = source.pointData.rt;

% Find indices of values greater than 1.5 and smaller than 0
sourceIndicesAbOutlier = find(sourceAb > 1.5 | sourceAb < 0);

% Find indices of values greater than 1.5 and smaller than 0
sourceIndicesRtOutlier = find(sourceRt > 1.5 | sourceRt < 0);

% Find the 4 nearest neighbors for all source points
sourceNeighborIndices = knnsearch(source.points, source.points, 'k', 4); % Include itself as the first neighbor

% Exclude the point itself (the first neighbor)
sourceNeighborIndices = sourceNeighborIndices(:, 2:end);

% Filter neighbors with values greater than 1.5
validNeighborsAb = sourceNeighborIndices(sourceAb(sourceNeighborIndices) < 1.5 & sourceAb(sourceNeighborIndices) >= 0);
% Filter neighbors with values greater than 1.5
validNeighborsRt = sourceNeighborIndices(sourceRt(sourceNeighborIndices) < 1.5 & sourceRt(sourceNeighborIndices) >= 0);


% Calculate the average values for valid neighbors
validNeighborValuesAb = sourceAb(validNeighborsAb);

% Calculate the average of the values for the 3 nearest neighbors if there are sourceIndices
if ~isempty(sourceIndicesAbOutlier)
    % Calculate the average values for sourceIndices
    sourceAveragesAb = zeros(size(sourceIndicesAbOutlier));
    sourceAveragesAb = accumarray(validNeighborsAb, validNeighborValuesAb, [], @mean);

    % Ensure that each sourceIndex has at least 3 valid neighbors
    invalidIndices = find(sourceAveragesAb(validNeighborsAb) < 1.5 & sourceAveragesAb(validNeighborsAb) >= 0);
    sourceIndicesAbOutlier(invalidIndices) = [];

    % Update sourceAb with the calculated averages
    sourceAb(sourceIndicesAbOutlier) = sourceAveragesAb(sourceIndicesAbOutlier);
end


% Filter neighbors with values greater than 1.5 or smaller than 0
validNeighborsRt = sourceNeighborIndices(sourceRt(sourceNeighborIndices) > 1.5 | sourceRt(sourceNeighborIndices) < 0);

% Calculate the average values for valid neighbors
validNeighborValuesRt = sourceRt(validNeighborsRt);

% Calculate the average of the values for the 3 nearest neighbors if there are sourceIndices
if ~isempty(sourceIndicesRtOutlier)
    % Calculate the average values for sourceIndices
    sourceAveragesRt = accumarray(validNeighborsRt, validNeighborValuesRt, [], @mean);

    % Ensure that each sourceIndex has at least 3 valid neighbors
    invalidIndices = find(sourceAveragesRt(validNeighborsRt) < 1.5 & sourceAveragesRt(validNeighborsRt) >= 0);
    sourceIndicesRtOutlier(invalidIndices) = [];

    % Update sourceRt with the calculated averages
    sourceRt(sourceIndicesRtOutlier) = sourceAveragesRt(sourceIndicesRtOutlier);
end

%keyboard;

% Update source and target data with modified ab values
source.pointData.ab = sourceAb;
source.pointData.rt = sourceRt;

%% recalculate sin and cos
source.pointData.rt = min(max(source.pointData.rt,0),1);
source.pointData.rtSin = sin(2*pi*source.pointData.rt);
source.pointData.rtCos = cos(2*pi*source.pointData.rt);

% repeat for target

% Process source.pointData.ab
targetAb = target.pointData.ab;
targetRt = target.pointData.rt;

% Find indices of values greater than 1.5 and smaller than 0
targetIndicesAbOutlier = find(targetAb > 1.5 | targetAb < 0);


% Find indices of values greater than 1.5 and smaller than 0
targetIndicesRtOutlier = find(targetRt > 1.5 | targetRt < 0);

% Find the 4 nearest neighbors for all target points
targetNeighborIndices = knnsearch(target.points, target.points, 'k', 4); % Include itself as the first neighbor

% Exclude the point itself (the first neighbor)
targetNeighborIndices = targetNeighborIndices(:, 2:end);



% Filter neighbors with values greater than 1.5
validNeighborsAb = targetNeighborIndices(targetAb(targetNeighborIndices) < 1.5 & targetAb(targetNeighborIndices) > 0);
% Filter neighbors with values greater than 1.5
validNeighborsRt = targetNeighborIndices(targetRt(targetNeighborIndices) < 1.5 & targetRt(targetNeighborIndices) > 0);


% Filter neighbors with values greater than 1.5 or smaller than 0
validNeighborsRt = targetNeighborIndices(targetRt(targetNeighborIndices) > 1.5 | targetRt(targetNeighborIndices) < 0);

% Calculate the average values for valid neighbors
validNeighborValuesRt = targetRt(validNeighborsRt);

% Calculate the average of the values for the 3 nearest neighbors if there are targetIndices
if ~isempty(targetIndicesRtOutlier)
    % Calculate the average values for targetIndices
    targetAveragesRt = accumarray(validNeighborsRt, validNeighborValuesRt, [], @mean);

    % Ensure that each targetIndex has at least 3 valid neighbors
    invalidIndices = find(targetAveragesRt(validNeighborsRt) < 1.5 & targetAveragesRt(validNeighborsRt) >= 0);
    targetIndicesRtOutlier(invalidIndices) = [];

    % Update targetRt with the calculated averages
    targetRt(targetIndicesRtOutlier) = targetAveragesRt(targetIndicesRtOutlier);
end


% Filter neighbors with values greater than 1.5 or smaller than 0
validNeighborsAb = targetNeighborIndices(targetAb(targetNeighborIndices) > 1.5 | targetAb(targetNeighborIndices) < 0);

% Calculate the average values for valid neighbors
validNeighborValuesAb = targetAb(validNeighborsAb);

% Calculate the average of the values for the 3 nearest neighbors if there are targetIndices
if ~isempty(targetIndicesAbOutlier)
    % Calculate the average values for targetIndices
    targetAveragesAb = accumarray(validNeighborsAb, validNeighborValuesAb, [], @mean);

    % Ensure that each targetIndex has at least 3 valid neighbors
    invalidIndices = find(targetAveragesAb(validNeighborsAb) < 1.5 & targetAveragesAb(validNeighborsAb) >= 0);
    targetIndicesAbOutlier(invalidIndices) = [];

    % Update targetAb with the calculated averages
    targetAb(targetIndicesAbOutlier) = targetAveragesAb(targetIndicesAbOutlier);
end


% Update target and target data with modified ab values
target.pointData.ab = targetAb;
target.pointData.rt = targetRt;

% based on modification, recalculate sin and cos for rt

target.pointData.rt = min(max(target.pointData.rt,0),1);
target.pointData.rtSin = sin(2*pi*target.pointData.rt);
target.pointData.rtCos = cos(2*pi*target.pointData.rt);



%% Scale ventricular coords to have a similar change across one tet.

if verbose, tic; fprintf('Scaling coordinates...                 '); end

tv_cells = round(source.pointData.tv(source.cells));
tv_norm = mean(vtkEdgeLengths(source)) / norm(max(source.points,[],1)-min(source.points,[],1));

rtSin_cells = source.pointData.rtSin(source.cells);
rtSin_norm = median(max(rtSin_cells,[],2) - min(rtSin_cells,[],2));

rtCos_cells = source.pointData.rtCos(source.cells);
rtCos_norm = median(max(rtCos_cells,[],2) - min(rtCos_cells,[],2));

ab_cells = source.pointData.ab(source.cells);
ab_norm = median(max(ab_cells,[],2) - min(ab_cells,[],2));

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

    X = NaN(size(source.cells,1), 5);
    X(:,1) = round(mean(tv_cells,2)) / tv_norm;
    X(:,2) = mean(ab_cells,2) / ab_norm;
    X(:,3) = mean(rtSin_cells,2) / rtSin_norm .* sqrt(mean(ab_cells,2));
    X(:,4) = mean(rtCos_cells,2) / rtCos_norm .* sqrt(mean(ab_cells,2));
    X(:,5) = mean(tm_cells,2) / tm_norm;


    Y = NaN(numTarPoints, 5);
    Y(:,1) = round(target.pointData.tv) / tv_norm;
    Y(:,2) = target.pointData.ab / ab_norm;
    Y(:,3) = target.pointData.rtSin / rtSin_norm .* sqrt(target.pointData.ab);
    Y(:,4) = target.pointData.rtCos / rtCos_norm .* sqrt(target.pointData.ab);
    Y(:,5) = target.pointData.tm / tm_norm;

    Mdl1 = KDTreeSearcher(X);

    pointIds = knnsearch(Mdl1, Y);

    if verbose, fprintf('%.1f seconds\n', toc); end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For each initial centroid, find all centroids within a radius
    % (euclidean distance in actual, cartesian coords).
    % This is done separately for left and right.

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
        source.pointData.rtSin .* sqrt(source.pointData.ab) ...
        source.pointData.rtCos .* sqrt(source.pointData.ab) ...
        source.pointData.tm ...
        ones(numSrcPoints,1) ...
        ]';

    t_coords = [ ...
        target.pointData.ab ...
        target.pointData.rtSin .* sqrt(target.pointData.ab) ...
        target.pointData.rtCos .* sqrt(target.pointData.ab) ...
        target.pointData.tm ...
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

    X = NaN(numSrcPoints, 5);
    X(:,1) = round(source.pointData.tv) / tv_norm;
    X(:,2) = source.pointData.ab / ab_norm;
    X(:,3) = source.pointData.rtSin / rtSin_norm .* sqrt(source.pointData.ab);
    X(:,4) = source.pointData.rtCos / rtCos_norm .* sqrt(source.pointData.ab);
    X(:,5) = source.pointData.tm / tm_norm;

    Y = NaN(numTarPoints, 5);
    Y(:,1) = round(target.pointData.tv) / tv_norm;
    Y(:,2) = target.pointData.ab / ab_norm;
    Y(:,3) = target.pointData.rtSin / rtSin_norm .* sqrt(target.pointData.ab);
    Y(:,4) = target.pointData.rtCos / rtCos_norm .* sqrt(target.pointData.ab);
    Y(:,5) = target.pointData.tm / tm_norm;

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

end
