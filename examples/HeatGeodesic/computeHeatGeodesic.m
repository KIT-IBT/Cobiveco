function [d,u] = computeHeatGeodesic(L, G, D, boundaryIds, t, tol, maxit)
% Computes an approx. solution phi of the eikonal equation |grad(phi)| = 1
% with boundary conditions phi(boundaryIds) = 0
% using the heat method (https://www.doi.org/10.1145/2516971.2516977).
%
% [d,u] = computeHeatGeodesic(L, G, D, boundaryIds, t, tol, maxit)
%
% Inputs:
%   L: Laplacian matrix computed with cotmatrix() of gptoolbox [numPoints x numPoints]
%   G: Gradient matrix computed with grad() of gptoolbox [3*numCells x numPoints]
%   D: Divergence matrix computed with div() of gptoolbox [numPoints x 3*numCells]
%   boundaryIds: IDs of starting points [numBoundaryPoints x 1]
%   t: time step used for solving the heat equation
%   tol: tolerance for bicgstab
%   maxit: maximum number of iterations for bicgstab
%
% Outputs:
%   d: Geodesic distance field [numPoints x 1]
%   u: Solution of the heat equation [numPoints x 1]
%
% Written by Steffen Schuler, Institute of Biomedical Engineering, KIT

N = length(L);

% set up initial conditions for heat equation
u0 = zeros(N,1);
u0(boundaryIds) = 1;

% solve heat equation
u = solveHeat(L, u0, t, tol, maxit);

% compute divergence of normalized gradient field of heat solution
X = -normalizedGradField(G, u, 1e-100, true);
div_X = 2*D*X(:);

% solve Poisson equation
boundaryVal = zeros(numel(boundaryIds),1);
d = solvePoisson(L, div_X, boundaryIds, boundaryVal, tol, maxit);

end