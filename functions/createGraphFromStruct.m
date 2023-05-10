function G = createGraphFromStruct(Struct)
% Extracts edges from a surface belonging to a matlab struct using "graph" function.
% Creates a graph based on edges and nodes in the surface using vtk functions.
% Automatically adds  a value for the node numbering in the table 
% allowing for consistent node numbering even when removing nodes or edges.
%
% G = createGraphFromStruct(Struct)
%
% Inputs:
%   Struct: tetrahedral volume mesh with nodes [numSourcePoints x 3] and elements/cells [numSourceCells x 4]
%
% Outputs:
%   G: Graph with table for node indices
%
% Written by Lisa Pankewitz


% extract Surface if needed

if isfield(Struct,'sur') && ~isempty(Struct.sur)
    fprintf('Surface already exists.\n')
else
    Struct.sur = vtkDataSetSurfaceFilter(vtkDeleteDataArrays(Struct.vol));
end

% extract edges using vtkExtractEdges

edgestruct = vtkExtractEdges(Struct.sur);

% add table value for consistent node numbering

node_indices = string(1:size(Struct.sur.points,1));

% compute edge lengths / distances to use as edge weights in the graph
vertex_a = edgestruct.cells(:,1);
vertex_b = edgestruct.cells(:,2);
p_a = edgestruct.points(vertex_a,:);
p_b = edgestruct.points(vertex_b,:);

norm_degree = 2; % L2 norm (aka. Eucledian distance)
axis = 2; % Compute norm over the xyz components
distances = vecnorm(p_a - p_b, norm_degree, axis);

% create graph with table for node indices
G = graph(vertex_a, vertex_b, distances, node_indices);


end