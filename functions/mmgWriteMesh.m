function mmgWriteMesh(mesh, filename)

fid = fopen(filename, 'w');
if ~fid
    error('Could not open %s for writing.', filename);
end

fprintf(fid, '%s\n', '# MESH #');
fprintf(fid, '%s\n\n', 'MeshVersionFormatted 2');
fprintf(fid, '%s\n\n', 'Dimension 3');

fprintf(fid, 'Vertices\n');
fprintf(fid, '%d\n', size(mesh.points,1));
if ~isfield(mesh, 'pointData') || ~isfield(mesh.pointData, 'class')
    mesh.pointData.class = zeros(size(mesh.points,1),1);
end
fprintf(fid, '%f %f %f %d\n', [mesh.points mesh.pointData.class]');
fprintf(fid, '\n');

fprintf(fid, 'Tetrahedra\n');
fprintf(fid, '%d\n', size(mesh.cells,1));
if ~isfield(mesh, 'cellData') || ~isfield(mesh.cellData, 'class')
    mesh.cellData.class = zeros(size(mesh.cells,1),1);
end
fprintf(fid, '%d %d %d %d\n\n', [mesh.cells mesh.cellData.class]');
fprintf(fid, '\n');

fprintf(fid, 'End \n');
fclose(fid);

end
