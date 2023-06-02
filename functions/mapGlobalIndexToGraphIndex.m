function [graphIndex] = mapGlobalIndexToGraphIndex(G, globalIndex, desc)
% Find node in graph corresponding to global index
%
% [graphIndex] = mapGlobalIndexToGraphIndex(G, globalIndex, desc)
%
% Input: 
%     G, Graph object
%     globalIndex, int
%     desc
%
% Output:
%     graphIndex, int
% 

% Node name is stored as string

globalIndexStr = string(globalIndex);
graphIndex = findnode(G, globalIndexStr);
if graphIndex == 0
    warning('Could not find global index %d (%s) in graph', globalIndex, desc);
end

end

