function M = ichol_autocomp(A, opts)

% Autocompiles the incomplete Cholesky factorization of A with 
% options specified by the structure options.
%
% M = ichol_autocomp(A, opts)
%
% Inputs:
%
%   A: square sparse matrix, double 
%   opts: options (see MATLAB help ichol)
%
% Output:
%
%   M: the factorized matrix, double 

opts.diagcomp = 1e-3*(max(sum(abs(A),2)./diag(A))-2);

for i = 1:50
    try
        M = ichol(A, opts);
        break;
    catch
        opts.diagcomp = 2*opts.diagcomp;
    end
end

if ~exist('M', 'var')
    error('Diagonal compensation for ichol failed.');
end

end