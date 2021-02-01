function splitMesh(o)

if ~o.available.transventricular
    o.computeTransventricular;
end
o.printStatus('Splitting mesh...                 ');
t = toc;

tmp = o.bv.vol;
tmp.pointData.ids = int32(1:size(tmp.points,1))';
tmp.pointData.class = repmat(uint8(0),size(tmp.points,1),1);
tmp.pointData.class(o.bv.surToVol) = o.bv.sur.pointData.class;
tmp.pointData.tv = o.bv.tv;

o.lv.vol = vtkTetConnectivityFilter(vtkThreshold(tmp, 'points', 'tv', [-inf 0.5]));
rid = double(mode(o.lv.vol.cellData.RegionId));
o.lv.vol = vtkThreshold(o.lv.vol, 'cells', 'RegionId', [rid rid]);
o.lv.ids = o.lv.vol.pointData.ids;
o.lv.sur = vtkDataSetSurfaceFilter(o.lv.vol);
o.lv.surToVol = vtkMapPointIds(o.lv.vol, o.lv.sur);
o.lv.vol = vtkDeleteDataArrays(o.lv.vol);
o.lv.sur = vtkDeleteDataArrays(o.lv.sur, {'pointData.class'});
o.lv.sur.cellData.class = min(o.lv.sur.pointData.class(o.lv.sur.cells),[],2);
o.lv.L = cotmatrix(double(o.lv.vol.points), double(o.lv.vol.cells));

o.rv.vol = vtkTetConnectivityFilter(vtkThreshold(tmp, 'points', 'tv', [0.5 inf]));
rid = double(mode(o.rv.vol.cellData.RegionId));
o.rv.vol = vtkThreshold(o.rv.vol, 'cells', 'RegionId', [rid rid]);
o.rv.ids = o.rv.vol.pointData.ids;
o.rv.sur = vtkDataSetSurfaceFilter(o.rv.vol);
o.rv.surToVol = vtkMapPointIds(o.rv.vol, o.rv.sur);
o.rv.vol = vtkDeleteDataArrays(o.rv.vol);
o.rv.sur = vtkDeleteDataArrays(o.rv.sur, {'pointData.class'});
o.rv.sur.cellData.class = min(o.rv.sur.pointData.class(o.rv.sur.cells),[],2);
o.rv.L = cotmatrix(double(o.rv.vol.points), double(o.rv.vol.cells));

if o.cfg.exportLevel > 2
    vtkWrite(o.lv.vol, sprintf('%slv.vtu', o.cfg.outPrefix));
    vtkWrite(o.rv.vol, sprintf('%srv.vtu', o.cfg.outPrefix));
end

o.printStatus(sprintf('%.1f seconds\n', toc-t), true);
o.available.splitMesh = true;

end