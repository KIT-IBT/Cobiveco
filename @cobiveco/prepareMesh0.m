function prepareMesh0(o)

o.printStatus('Preparing mesh 0...');
t = toc;

o.m0.vol = vtkRead(sprintf('%s.vtu', o.cfg.inPrefix));
o.m0.vol = vtkDeleteDataArrays(o.m0.vol);
o.m0.sur = vtkDataSetSurfaceFilter(o.m0.vol);       % surface of volume mesh
o.m0.surToVol = vtkMapPointIds(o.m0.vol, o.m0.sur); % mapping from surface to volume point ids
o.m0.meanEdgLen = mean(vtkEdgeLengths(o.m0.vol));   % mean edge length
o.result = o.m0.vol;

vtpFile = sprintf('%s.vtp', o.cfg.inPrefix);
if exist(vtpFile, 'file')
    % map classes provided in vtp file to surface of volume mesh
    sur = vtkRead(vtpFile);
    sur.pointData.class(sur.pointData.class>4) = 2;
    ids = vtkMapPointIds(sur, o.m0.sur);
    o.m0.sur.pointData.class = sur.pointData.class(ids);
else
    % map given boundary surfaces to surface of volume mesh and annotate classes
    surBase = vtkRead(sprintf('%s_base.ply', o.cfg.inPrefix));
    surEpi = vtkRead(sprintf('%s_epi.ply', o.cfg.inPrefix));
    surEndoLv = vtkRead(sprintf('%s_endo_lv.ply', o.cfg.inPrefix));
    surEndoRv = vtkRead(sprintf('%s_endo_rv.ply', o.cfg.inPrefix));
    surAppended = vtkAppendPolyData({surBase, surEpi, surEndoLv, surEndoRv});
    cnp = cumsum([size(surBase.points,1); size(surEpi.points,1); size(surEndoLv.points,1)]);
    ids = vtkMapPointIds(surAppended, o.m0.sur);
    o.m0.sur.pointData.class = repmat(uint8(0),size(o.m0.sur.points,1),1);
    o.m0.sur.pointData.class(ids<=cnp(1)) = 1;               % base
    o.m0.sur.pointData.class(ids>cnp(1) & ids<=cnp(2)) = 2;  % epi
    o.m0.sur.pointData.class(ids>cnp(2) & ids<=cnp(3)) = 3;  % endoLv
    o.m0.sur.pointData.class(ids>cnp(3)) = 4;                % endoRv
    mappingDist = sqrt(sum((o.m0.sur.points-surAppended.points(ids,:)).^2,2));
    o.m0.sur.pointData.class(mappingDist > o.cfg.mappingTol*o.m0.meanEdgLen) = 0;
end

P0 = double(o.m0.vol.points);
C0 = double(o.m0.vol.cells);
o.m0.L = cotmatrix(P0, C0);

if o.cfg.exportLevel > 1
    o.m0.debug = o.m0.vol;
    o.m0.debug.pointData.surClass = repmat(uint8(0),size(o.m0.debug.points,1),1);
    o.m0.debug.pointData.surClass(o.m0.surToVol) = o.m0.sur.pointData.class;
end

o.printStatus(sprintf('%.1f seconds\n', toc-t), true);
o.available.mesh0 = true;

end