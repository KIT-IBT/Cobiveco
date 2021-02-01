function u = solveHeat(L, u0, timestep, tol, maxit)
% Computes the solution u(t=timestep) of the heat equation L*u = du/dt
% with initial conditions u(t=0) = u0.
% Time integration is done using a single backward Euler step, resulting in
% (I-timestep*L)*u(t=timestep) = u0.
%
% u = solveHeat(L, u0, timestep, tol, maxit)
%
% Inputs:
%   L: Laplacian matrix computed with cotmatrix() of gptoolbox [numPoints x numPoints]
%   u0: Initial temperature distribution [numPoints x 1]
%   timestep: time step used for one single backward Euler step
%   tol: tolerance for bicgstab
%   maxit: maximum number of iterations for bicgstab
%
% Outputs:
%   u: Solution of the heat equation [numPoints x 1]
%
% Written by Steffen Schuler, Institute of Biomedical Engineering, KIT

if nargin < 6
    maxit = 500;
end
if nargin < 5
    tol = [];
end

N = length(L);

% set up coefficient matrix
A = speye(N)-timestep*L;

% apply reverse Cuthill-McKee reordering
p = symrcm(A);
A = A(p,p);
u0 = u0(p);

% compute incomplete LU factorization
[lowTri, upTri] = ilu(A);

% solve linear system
[x, flag, relres, iter] = bicgstab(A, u0, tol, maxit, lowTri, upTri);
if flag
    error('bicgstab failed at iteration %i with flag %i and relative residual %.1e.', iter, flag, relres);
end
u = NaN(N,1);
u(p) = x;

end