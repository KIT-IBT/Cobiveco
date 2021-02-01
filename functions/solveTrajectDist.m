function d = solveTrajectDist(G, T, boundaryIds, boundaryVal, tol, maxit)
% Computes the solution d of the "trajectory distance" equation dot(G*d,T) = 1
% with boundary conditions u(boundaryIds) = boundaryVal
% (see https://doi.org/10.1109/TMI.2003.817775, equation 3, where L0 = d).
%
% d = solveTrajectDist(G, T, boundaryIds, boundaryVal, tol, maxit)
%
% Inputs:
%   G: gradient operator matrix computed with grad() of gptoolbox [3*numCells x numPoints]
%   T: normalized tangent field defining the course of trajectories [numCells x 3]
%   boundaryIds: point IDs for Dirichlet boundary conditions [numBoundaryPoints x 1]
%   boundaryVal: values for Dirichlet boundary conditions [numBoundaryPoints x 1]
%   tol: tolerance for minres
%   maxit: maximum number of iterations for minres
%
% Outputs:
%   d: distance along trajectories (starting from boundaryIds with boundaryVal) [numPoints x 1]
%
% Written by Steffen Schuler, Institute of Biomedical Engineering, KIT

if nargin < 6
    maxit = 2000;
end
if nargin < 5
    tol = [];
end

if numel(boundaryIds) ~= numel(unique(boundaryIds))
    error('Duplicate entries found in boundaryIds.');
end

[nC,nDim] = size(T);

% For each cell, T_mat computes the dot product of the tangent field T
% with another vector field defined at the cells and represented as a
% 3*nC x 1 vector [x1; x2; ...; y1; y2; ...; z1; z2; ...]
i = reshape(repmat(1:nC, nDim, 1), [], 1);
j = reshape(reshape(1:nDim*nC, [], 3)', [], 1);
v = reshape(T', [], 1);
T_mat = sparse(i, j, v, nC, nDim*nC);

% For each cell, S computes the dot product of the tangent field T
% with the gradient of a scalar field defined at the points,
% i.e. a nP x 1 vector
S = T_mat * G;

% As nC > nP, the linear system S*d = 1 is overdetermined,
% so we minimize norm(S*d-1) by solving A*d = b
A = S'*S;
b = S'*ones(nC,1);

N = length(A);
K = numel(boundaryIds);
boundaryIds = double(boundaryIds);

% set up right-hand side
b = b - A(:,boundaryIds)*boundaryVal(:);
b(boundaryIds) = boundaryVal;

% add boundary conditions to coefficient matrix
A(boundaryIds,:) = sparse(1:K, boundaryIds, ones(K,1), K, N);
A(:,boundaryIds) = sparse(boundaryIds, 1:K, ones(K,1), N, K);

% apply reverse Cuthill-McKee reordering
p = symrcm(A);
A = A(p,p);
b = b(p);

icMat = ichol_autocomp(A);
[x, flag, relres, iter] = minres(A, b, tol, maxit, icMat, icMat');
if flag
    warning('minres failed at iteration %i with flag %i and relative residual %.1e.', iter, flag, relres);
end
d = NaN(N,1);
d(p) = x;

end