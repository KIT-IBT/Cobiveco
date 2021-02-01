function [mesh,status,cmdout] = mmg(mesh, sol, paramString)

[tmpdir,name] = fileparts(tempname);
meshfile = sprintf('%s/%s.mesh', tmpdir, name);
solfile = sprintf('%s/%s.sol', tmpdir, name);

mmgWriteMesh(mesh, meshfile);
mmgWriteSol(sol, solfile);

mpath = fileparts(mfilename('fullpath'));
[status,cmdout] = system(sprintf('%s/../dependencies/mmg/build/bin/mmg3d_O3 %s %s -sol %s %s', mpath, meshfile, meshfile, solfile, paramString));

mesh = struct();
if status==0
    mesh = mmgReadMesh(meshfile);
end

end
