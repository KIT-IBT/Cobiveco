function computeTransventricular(o)

% Computes the transventricular coordinate of the elements in the object, 
% i.e. it indicates if the element belongs to the left (0) or right (1) ventricle.
% The coordinate value is added in the field "tv" to the output objects.
% 
% computeTransventricular(o)
%
% Input: 
%   o: instance of the class cobiveco
% 
% Output: 
%   o.m0, object of class cobiveco [struct] (For details see cobiveco class documentation)


if ~o.available.mesh0
    o.prepareMesh0;
end
o.printStatus('Computing transventricular coordinate...');
t = toc;

idsEndoLv = o.m0.surToVol(o.m0.sur.pointData.class==3);
idsEndoRv = o.m0.surToVol(o.m0.sur.pointData.class==4);

ids = [idsEndoLv; idsEndoRv];
val = [zeros(size(idsEndoLv)); ones(size(idsEndoRv))];
o.m0.tvLaplace = solveLaplace(o.m0.L, ids, val, o.cfg.tol, o.cfg.maxit);
o.m0.tv = round(o.m0.tvLaplace);
o.result.pointData.tv = single(o.m0.tv);

if o.cfg.exportLevel > 1
    o.m0.debug.pointData.tvLaplace = single(o.m0.tvLaplace);
    o.m0.debug.pointData.tv = single(o.m0.tv);
    vtkWrite(o.m0.debug, sprintf('%sdebug0.vtu', o.cfg.outPrefix));
end

o.printStatus(sprintf('%.1f seconds\n', toc-t), true);
o.available.transventricular = true;

end