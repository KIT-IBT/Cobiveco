function M = ichol_autocomp(A, opts)

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