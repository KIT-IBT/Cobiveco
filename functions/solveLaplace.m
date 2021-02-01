function u = solveLaplace(L, boundaryIds, boundaryVal, tol, maxit)
% Computes the solution u of the Laplace equation L*u = 0
% with boundary conditions u(boundaryIds) = boundaryVal.
%
% u = solveLaplace(L, boundaryIds, boundaryVal, tol, maxit)
%
% Inputs:
%   L: Laplacian matrix computed with cotmatrix() of gptoolbox [numPoints x numPoints]
%   boundaryIds: point IDs for Dirichlet boundary conditions [numBoundaryPoints x 1]
%   boundaryVal: values for Dirichlet boundary conditions [numBoundaryPoints x 1]
%   tol: tolerance for bicgstab
%   maxit: maximum number of iterations for bicgstab
%
% Outputs:
%   u: Laplace solution [numPoints x 1]
%
% Written by Steffen Schuler, Institute of Biomedical Engineering, KIT

if nargin < 5
    maxit = 500;
end
if nargin < 4
    tol = [];
end

if numel(boundaryIds) ~= numel(unique(boundaryIds))
    error('Duplicate entries found in boundaryIds.');
end

N = length(L);
K = numel(boundaryIds);
boundaryIds = double(boundaryIds);

% set up right-hand side
b = -L(:,boundaryIds)*boundaryVal(:);
b(boundaryIds) = boundaryVal;

% add boundary conditions to coefficient matrix
L(boundaryIds,:) = sparse(1:K, boundaryIds, ones(K,1), K, N);
L(:,boundaryIds) = sparse(boundaryIds, 1:K, ones(K,1), N, K);

% apply reverse Cuthill-McKee reordering
p = symrcm(L);
L = L(p,p);
b = b(p);

% compute incomplete LU factorization
[lowTri, upTri] = ilu(L);

% solve linear system
[x, flag, relres, iter] = bicgstab(L, b, tol, maxit, lowTri, upTri);
if flag
    warning('bicgstab failed at iteration %i with flag %i and relative residual %.1e.', iter, flag, relres);
end
u = NaN(N,1);
u(p) = x;

end