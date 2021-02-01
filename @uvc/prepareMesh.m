function prepareMesh(o)

o.printStatus('Preparing mesh...');
t = toc;

o.bv.vol = vtkRead(strrep(sprintf('%s.vtu', o.cfg.inPrefix), '_for_UVC', ''));
o.bv.vol = vtkDeleteDataArrays(o.bv.vol);
o.bv.sur = vtkDataSetSurfaceFilter(o.bv.vol);       % surface of volume mesh
o.bv.surToVol = vtkMapPointIds(o.bv.vol, o.bv.sur); % mapping from surface to volume point ids
o.meanEdgLen = mean(vtkEdgeLengths(o.bv.vol));      % mean edge length

vtpFile = sprintf('%s.vtp', o.cfg.inPrefix);
sur = vtkRead(vtpFile);
ids = vtkMapPointIds(sur, o.bv.sur);
o.bv.sur.pointData.class = sur.pointData.class(ids);

o.bv.L = cotmatrix(double(o.bv.vol.points), double(o.bv.vol.cells));

o.result = o.bv.vol;

o.printStatus(sprintf('%.1f seconds\n', toc-t), true);
o.available.mesh = true;

end