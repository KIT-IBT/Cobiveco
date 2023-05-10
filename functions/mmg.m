function [mesh,status,cmdout] = mmg(mesh, sol, paramString)

% [mesh,status,cmdout] = mmg(mesh, sol, paramString)
% 
% Uses mmg meshing Software 
% (for more infos, see https://www.mmgtools.org/)


[tmpdir,name] = fileparts(tempname);
meshfile = sprintf('%s/%s.mesh', tmpdir, name);
solfile = sprintf('%s/%s.sol', tmpdir, name);

mmgWriteMesh(mesh, meshfile);
mmgWriteSol(sol, solfile);

mpath = fileparts(mfilename('fullpath'));
mmg_executable_path = sprintf('%s/../dependencies/mmg/build/bin/mmg3d_O3', mpath);
[status,cmdout] = system(sprintf('"%s" %s %s -sol %s %s', mmg_executable_path, meshfile, meshfile, solfile, paramString));
mesh = struct();
if status==0
    mesh = mmgReadMesh(meshfile);
elseif status == 126 || status == 127
    message = sprintf('failed to run mmg3d (file at "%s")', mmg_executable_path);
    error(message);
    error(cmdout);
end

end
