function computeTransmural(o)

if ~o.available.splitMesh
    o.splitMesh;
end
o.printStatus('Computing transmural coordinate...');
t = toc;

idsEndoLv = o.lv.surToVol(o.lv.sur.pointData.class==3);
idsEpiLv = o.lv.surToVol(o.lv.sur.pointData.class~=1 & o.lv.sur.pointData.class~=3);
ids = [idsEndoLv; idsEpiLv];
val = [zeros(size(idsEndoLv)); ones(size(idsEpiLv))];
o.lv.tm = solveLaplace(o.lv.L, ids, val, o.cfg.tol, o.cfg.maxit);

idsEndoRv = o.rv.surToVol(o.rv.sur.pointData.class==4);
idsEpiRv = o.rv.surToVol(o.rv.sur.pointData.class==2);
ids = [idsEndoRv; idsEpiRv];
val = [zeros(size(idsEndoRv)); ones(size(idsEpiRv))];
o.rv.tm = solveLaplace(o.rv.L, ids, val, o.cfg.tol, o.cfg.maxit);

o.bv.tm = NaN(size(o.bv.vol.points,1),1);
o.bv.tm(o.lv.ids) = o.lv.tm;
o.bv.tm(o.rv.ids) = o.rv.tm;
bvNan = isnan(o.bv.tm);
lvNan = bvNan & o.bv.tv==0;
rvNan = bvNan & o.bv.tv==1;
lvNotNan = find(~bvNan & o.bv.tv==0);
rvNotNan = find(~bvNan & o.bv.tv==1);
o.bv.tm(lvNan) = o.bv.tm(lvNotNan(knnsearch(o.bv.vol.points(lvNotNan,:), o.bv.vol.points(lvNan,:))));
o.bv.tm(rvNan) = o.bv.tm(rvNotNan(knnsearch(o.bv.vol.points(rvNotNan,:), o.bv.vol.points(rvNan,:))));

o.result.pointData.tm = single(o.bv.tm);

if o.cfg.exportLevel > 2
    vtkWrite(o.result, sprintf('%sresult.vtu', o.cfg.outPrefix));
end

o.printStatus(sprintf('%.1f seconds\n', toc-t), true);
o.available.transmural = true;

end