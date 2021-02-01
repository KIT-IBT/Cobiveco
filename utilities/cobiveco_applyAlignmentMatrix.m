function vtk = cobiveco_applyAlignmentMatrix(vtk, R)
% Transforms a heart mesh using an alignment matrix as computed by cobiveco.
% After alignment, the heart axes should be parallel to x,y,z and the heart
% should be centered at (0,0,0).
%
% Syntax:
%  vtk = cobiveco_applyAlignmentMatrix(vtk, R)
%
% Inputs:
% - vtk: VTK struct of the heart mesh
% - R:   alignment matrix as computed by cobiveco
%
% Output:
% - vtk: VTK struct of the transformed heart mesh
%
% Written in 2020 by Steffen Schuler
% Institute of Biomedical Engineering, KIT
% www.ibt.kit.edu

vtk.points = [vtk.points ones(size(vtk.points,1),1)]*R';
vtk.points(:,end) = [];

end