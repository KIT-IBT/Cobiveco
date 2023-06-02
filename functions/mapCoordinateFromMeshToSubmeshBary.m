

function [struct2] = mapCoordinateFromMeshToSubmeshBary(struct1, struct2, coord1, coord2)
    % Using Barycentric Interpolation to map coordinates from one structe to another, a substructure of a mesh.
    %
    % [struct2] = mapCoordinateFromMeshToSubmeshBary(struct1, struct2, coord1, coord2)
    %
    % Inputs:
    %
    %   struct1 [numPoints x 3]
    %   struct2 [numPoints x 3]
    %   coord1, str with value "tm" or "tv"; double [numPointsx1], except for tv  int [numPointsx1]
    %   coord2, str with value "tm" or "tv"; double [numPointsx1], except for tv  int [numPointsx1]
    %
    % Output:
    %
    %   struct2 [numPoints x 3]
    %
    % Written by Lisa Pankewitz
    
    if nargin < 4
        if coord1 == 'tv'
            Matrix = baryInterpMat(struct1.vol.points,struct1.vol.cells, struct2.vol.points);
            struct2.vol.pointData.tv = Matrix * struct1.tv;
        elseif coord1 == 'tm'
            Matrix = baryInterpMat(struct1.vol.points,struct1.vol.cells, struct2.vol.points);
            struct2.vol.pointData.tm = Matrix * struct1.tm;

        else
            disp("Please enter a valid coordinate. Aborting.")
        end
        
    else
        if coord1 == 'tv' & coord2 == 'tm'
            Matrix = baryInterpMat(struct1.vol.points,struct1.vol.cells, struct2.vol.points);
            struct2.vol.pointData.tv = Matrix * struct1.tv;
            struct2.vol.pointData.tm = Matrix * struct1.tm;


        elseif coord1 == 'tm' & coord2 == 'tv'
            Matrix = baryInterpMat(struct1.vol.points,struct1.vol.cells, struct2.vol.points);
            struct2.vol.pointData.tv = Matrix * struct1.tv;
            struct2.vol.pointData.tm = Matrix * struct1.tm;

        else
            disp("Please enter a valid coordinate option. You can enter either 'tv' or 'tm'. Aborting.")
        end

    
    end
    
    end
    

