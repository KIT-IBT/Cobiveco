function computeTransventricular(o)

if ~o.available.mesh
    o.prepareMesh;
end
o.printStatus('Computing transventricular coordinate...');
t = toc;

idsEndoLv = o.bv.surToVol(o.bv.sur.pointData.class==3);
idsEndoRv = o.bv.surToVol(o.bv.sur.pointData.class==4);
idsSeptumRv = o.bv.surToVol(o.bv.sur.pointData.class==6);

ids = [idsEndoLv; idsSeptumRv; idsEndoRv];
val = [zeros(size(idsEndoLv)); repmat(0.5, size(idsSeptumRv)); ones(size(idsEndoRv))];
o.bv.tv = round(solveLaplace(o.bv.L, ids, val, o.cfg.tol, o.cfg.maxit)-0.01);

o.result.pointData.tv = single(o.bv.tv);

if o.cfg.exportLevel > 2
    vtkWrite(o.result, sprintf('%sresult.vtu', o.cfg.outPrefix));
end

o.printStatus(sprintf('%.1f seconds\n', toc-t), true);
o.available.transventricular = true;

end