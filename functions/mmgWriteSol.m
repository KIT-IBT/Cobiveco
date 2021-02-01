function mmgWriteSol(sol, filename)

fid = fopen(filename, 'w');
if ~fid
    error('Could not open %s for writing.', filename);
end

fprintf(fid, '%s\n', '# SOL #');
fprintf(fid, '%s\n\n', 'MeshVersionFormatted 2');
fprintf(fid, '%s\n\n', 'Dimension 3');

fprintf(fid, 'SolAtVertices\n');
fprintf(fid, '%d\n', size(sol,1));

fprintf(fid, '1 1\n');
fprintf(fid, '%1.15e\n', sol');
fprintf(fid, '\n');

fprintf(fid, 'End \n');
fclose(fid);

end
