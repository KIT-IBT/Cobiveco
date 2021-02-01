function T = normalizedGradField(G, u, thresh, warn, P, C)
% Computes the normalized gradient field T of a scalar field u.
% Gradient vectors in T are normalized to length 1.
%
% T = normalizedGradField(G, u, thresh, warn, P, C)
%
% Inputs:
%   G: gradient operator matrix computed with grad() of gptoolbox [3*numCells x numPoints]
%   u: scalar field [numPoints x 1]
%   thresh: threshold on the gradient length to avoid division by zero during normalization
%   warn: whether to print a warning when the gradient is too small
%   P, C: point list [numPoints x 3] and cell list [numCells x 4]
%         If these lists are provided, the nearest neighbor gradients are
%         used for cells with gradient lengths below the threshold.
%         Otherwise, the normalization is simply skipped for these cells.
%
% Outputs:
%   T: Normalized gradient field [numCells x 3]
%
% Written by Steffen Schuler, Institute of Biomedical Engineering, KIT

if nargin < 4 || isempty(warn)
    warn = true;
end

T = reshape(G*u, [], 3);
len = sqrt(sum(T.^2,2));
belowThresh = len < thresh;
len(belowThresh) = 1;
T = T./len;

if any(belowThresh)
    if nargin > 5
        if warn
            warning('Gradient length is smaller than threshold in the following cells - using gradient of nearest neighbor above threshold for these cells.\n%s', sprintf('%i ', find(belowThresh)));
        end
        centroids = squeeze(mean(reshape(P(C,:),[],size(C,2),size(P,2)),2));
        ids = find(~belowThresh);
        ids = ids(knnsearch(centroids(~belowThresh,:), centroids(belowThresh,:)));
        T(belowThresh,:) = T(ids,:);
    else
        if warn
            warning('Gradient length is smaller than threshold in the following cells - skipped normalization for these cells.\n%s', sprintf('%i ', find(belowThresh)));
        end
    end
end

end