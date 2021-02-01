function mesh = mmgReadMesh(filename)

fid = fopen(filename, 'r');
if ~fid
    error('Could not open %s for reading.', filename);
end

str = fgets(fid);
while str ~= -1
    
    if startsWith(str,'Vertices')
        N = str2double(fgets(fid));
        points = fscanf(fid,'%f %f %f %d',[4 N])';
        
    elseif startsWith(str,'Tetrahedra')
        N = str2double(fgets(fid));
        cells = fscanf(fid,'%d %d %d %d %d',[5 N])';
        
    end
    
    str = fgets(fid);
end

fclose(fid);

mesh.points = points(:,1:3);
mesh.cells = int32(cells(:,1:4));
mesh.cellData.class = uint8(cells(:,5));
mesh.cellTypes = repmat(uint8(10), size(cells,1), 1);

end
