function surnormal = pcaSurfaceNormal(sur, ref)
    % Calculate the surface normal of the input struct using the third principal component.
    %
    % surnormal = pcaSurfaceNormal(sur, ref)
    %
    % Inputs:
    %
    %   sur struct, surface representing the valve annuli with points being [numSurPoints x 3]
    %   ref, reference vector with known direction, double [1x3]
    %
    % Outputs:
    %
    %   surnormal, double [1x3]
    %
    % Written by Lisa Pankewitz

    centroidsSur = squeeze(mean(reshape(sur.points(sur.cells,:),[],size(sur.cells,2),size(sur.points,2)),2));
    areaSur = doublearea(sur.points, sur.cells)/2;
    pc = pca(centroidsSur, 'Weights', areaSur./mean(areaSur));
    surnormal = pc(:,end)';
    
    % check the dot product between the reference and the plane normal 1, if negative, flip
    if  dot(ref,surnormal) < 0
        surnormal = -surnormal;
    end



end
    
    
    
    