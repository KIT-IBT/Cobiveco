function target = updateOutliersWithNeighborValues(target)
    % Process target.pointData.Ab
    targetAb = target.pointData.ab;

    % Find indices of values greater than 1.5 and smaller than 0
    targetIndicesAbOutlier = find(targetAb > 1.5 | targetAb < 0);


    % Find the 4 nearest neighbors for all target points
    targetNeighborIndices = knnsearch(target.points, target.points(targetIndicesAbOutlier, :), 'k', 9); % Include itself as the first neighbor

    % Exclude the point itself (the first neighbor)
    current_point_index = targetNeighborIndices(:, 1);
    targetNeighborIndices = targetNeighborIndices(:, 2:end);

    % Calculate the average of the values for the 3 nearest neighbors if there are targetIndices
    if ~isempty(targetIndicesAbOutlier)
        % Calculate the average values for targetIndices
        targetAveragesAb = zeros(size(targetIndicesAbOutlier));

        % Iterate over each row in targetNeighborIndices
        for i = 1:size(targetNeighborIndices, 1)
            % Get the indices of the 4 nearest neighbors
            neighbors = targetNeighborIndices(i, :);

            % Initialize an array to store valid neighbors
            validNeighbors = [];

            % Check the criteria for the first three values
            for j = 1:8
                if targetAb(neighbors(j)) <= 1.5 && targetAb(neighbors(j)) >= 0
                    % If the value meets the criteria, add it to validNeighbors
                    validNeighbors = [validNeighbors, targetAb(neighbors(j))];
                end
            end

            % Check if there are at least 3 valid neighbors
            if length(validNeighbors) >= 3
                % Calculate the average based on the first, second, and third valid neighbors
                targetAveragesAb(i) = mean(validNeighbors([1, 2, 3]));

            elseif length(validNeighbors) >= 2
                % Calculate the average based on the first, second, and third valid neighbors
                targetAveragesAb(i) = mean(validNeighbors([1, 2]));
                disp(targetAveragesAb(i))

            else
                disp(targetAveragesAb(current_point_index(i)));
                targetAveragesAb(i) = validNeighbors(1);
                disp("Cannot Determine at least 2 neighbors. Check manually.")
            end
        end
    end

    if ~isempty(targetIndicesAbOutlier)
        % Update targetAb with the calculated averages
        targetAb(targetIndicesAbOutlier) = targetAveragesAb;


    end

    % Update the target structure
    target.pointData.ab = targetAb;

    % Process target.pointData.Rt
    targetRt = target.pointData.rt;

    % Find indices of values greater than 1.5 and smaller than 0
    targetIndicesRtOutlier = find(targetRt > 1.5 | targetRt < 0);

    % Find the 4 nearest neighbors for all target points
    targetNeighborIndices = knnsearch(target.points, target.points(targetIndicesRtOutlier, :), 'k', 9); % Include itself as the first neighbor

    % Exclude the point itself (the first neighbor)
    current_point_index = targetNeighborIndices(:, 1);
    targetNeighborIndices = targetNeighborIndices(:, 2:end);

    % Calculate the average of the values for the 3 nearest neighbors if there are targetIndices
    if ~isempty(targetIndicesRtOutlier)
        % Calculate the average values for targetIndices
        targetAveragesRt = zeros(size(targetIndicesRtOutlier));

        % Iterate over each row in targetNeighborIndices
        for i = 1:size(targetNeighborIndices, 1)
            % Get the indices of the 4 nearest neighbors
            neighbors = targetNeighborIndices(i, :);

            % Initialize an array to store valid neighbors
            validNeighbors = [];

            % Check the criteria for the first three values
            for j = 1:8
                if targetRt(neighbors(j)) <= 1.5 && targetRt(neighbors(j)) >= 0
                    % If the value meets the criteria, add it to validNeighbors
                    validNeighbors = [validNeighbors, targetRt(neighbors(j))];
                    disp(validNeighbors);
                    disp(length(validNeighbors))
                end
            end

            % Check if there are at least 3 valid neighbors
            if length(validNeighbors) >= 3
                % Calculate the average based on the first, second, and third valid neighbors
                targetAveragesRt(i) = mean(validNeighbors([1, 2, 3]));
                disp(targetAveragesRt(i))

            elseif length(validNeighbors) >= 2
                % Calculate the average based on the first, second, and third valid neighbors
                targetAveragesRt(i) = mean(validNeighbors([1, 2]));
                disp(targetAveragesRt(i))

            else
                disp(targetRt(current_point_index(i)));
                targetAveragesRt(i) = validNeighbors(1);
                disp("Cannot Determine at least 2 neighbors. Check manually.")

            end
        end
    end

    if ~isempty(targetIndicesRtOutlier)
        % Update targetRt with the calculated averages
        targetRt(targetIndicesRtOutlier) = targetAveragesRt;
    end

    % Update the target structure
    target.pointData.rt = targetRt;

end