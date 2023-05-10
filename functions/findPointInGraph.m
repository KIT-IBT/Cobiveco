
function [idxInGraphPoint1, idxInGraphPoint2] = findPointInGraph(struct, point1, point2, excludedPoints)

% Finds the point corresponding to the node in the graph
%
% [idxInGraphPoint1, idxInGraphPoint2] = findPointInGraph(struct, point1, point2, excludedPoints)
%
% Inputs:
%
%   struct [numSurPoints x 3]
%   point1, double [1x3]
%   point2, double [1x3]
%   excludedPoints
%
% Output:
%
%    idxInGraphPoint1, int [1x1]
%    idxInGraphPoint2, int [1x1]
%
% Written by Lisa Pankewitz

if nargin < 4
    point1InGraph = ismembertol(struct.points,point1,0.1,'ByRows',true,'DataScale',1.0);
    [idxInGraphPoint1,values] = find(point1InGraph==1);

    point2InGraph = ismembertol(struct.points,point2,0.1,'ByRows',true,'DataScale',1.0);
    [idxInGraphPoint2,values] = find(point2InGraph==1);

else
    tol = 0.1;
    numTries = 7;
    for i = 1:numTries
        point1InGraph = ismembertol(struct.points,point1,tol,'ByRows',true,'DataScale',1.0);
        % prevent matching any excluded points
        point1InGraph(excludedPoints) = 0;
        [idxInGraphPoint1,values] = find(point1InGraph==1);

        if ~isempty(idxInGraphPoint1)
            fprintf('Point 1 found with tol %g\n', tol);
            break;

        elseif i < numTries
            warning('Finding point in graph with tolerance of %.3f failed . Trying again with a slightly larger tolerance.', tol);
            tol = tol + 0.2;

        else
            error('identifying point failed ultimately after trying %i different thresholds (system command status %i).', numTries);

        end

    end

    for i = 1:numTries
        point2InGraph = ismembertol(struct.points,point2,tol,'ByRows',true,'DataScale',1.0);
        % prevent matching any excluded points
        point2InGraph(excludedPoints) = 0;
        [idxInGraphPoint2,values] = find(point2InGraph==1);

        if ~isempty(idxInGraphPoint2)
            fprintf('Point 2 found with tol %g\n', tol);
            break;
        elseif i < numTries
            warning('Finding point in graph with tolerance of %.3f failed . Trying again with a slightly larger tolerance.', tol);
            tol = tol + 0.2;

        else
            error('identifying point failed ultimately after trying %i different thresholds (system command status %i).', numTries);

        end

    end

    idxInGraphPoint1 = idxInGraphPoint1(1);
    idxInGraphPoint2 = idxInGraphPoint2(1);



end

end
