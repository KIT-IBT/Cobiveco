function prepareMesh0(o)

% Maps given boundary surfaces to surface of volume mesh and annotates classes.
% For the given object, the tetrahedral volume and the triangular boundary surfaces 
% are coregistered in the structure m0.
% The values are added in the fields "vol", "sur", "surToVol", "meanEdgLen" and "L" of the output object.
%
% prepareMesh0(o)
%
% Input: 
%   o: instance of the class cobiveco
% 
% Outputs: 
%   o.m0, object of class cobiveco [struct] (For details see cobiveco class documentation)


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
    epiBaseFile = sprintf('%s_epi_base.ply', o.cfg.inPrefix);
    if exist(epiBaseFile, 'file')
        surEpiBase = vtkRead(sprintf('%s_epi_base.ply', o.cfg.inPrefix));
        surEpiNoBase= vtkRead(sprintf('%s_epi_no_base.ply', o.cfg.inPrefix));
        o.surMv = vtkRead(sprintf('%s_mv.ply', o.cfg.inPrefix));
        o.surTv = vtkRead(sprintf('%s_tv.ply', o.cfg.inPrefix));
        o.surAv = vtkRead(sprintf('%s_av.ply', o.cfg.inPrefix));
        o.surPv = vtkRead(sprintf('%s_pv.ply', o.cfg.inPrefix));
        surEndoLv = vtkRead(sprintf('%s_endo_lv.ply', o.cfg.inPrefix));
        surEndoRv = vtkRead(sprintf('%s_endo_rv.ply', o.cfg.inPrefix));
        surAppended = vtkAppendPolyData({surEpiBase, surEpiNoBase, surEndoLv, surEndoRv, o.surMv,o.surTv, o.surAv, o.surPv});
        cnp = cumsum([size(surEpiBase.points,1); size(surEpiNoBase.points,1);size(surEndoLv.points,1); size(surEndoRv.points,1);size(o.surMv.points,1); size(o.surTv.points,1); size(o.surAv.points,1); size(o.surPv.points,1)]);
        ids = vtkMapPointIds(surAppended, o.m0.sur);
        o.m0.sur.pointData.class = repmat(uint8(0),size(o.m0.sur.points,1),1);    
        o.m0.sur.pointData.class(ids<=cnp(1)) = 1;               % EpiBase
        o.m0.sur.pointData.class(ids>cnp(1) & ids<=cnp(2)) = 2;  % surEpiNoBase
        o.m0.sur.pointData.class(ids>cnp(2) & ids<=cnp(3)) = 3;  % LV
        o.m0.sur.pointData.class(ids>cnp(3) & ids<=cnp(4)) = 4;  % Rv
        o.m0.sur.pointData.class(ids>cnp(4) & ids<=cnp(5)) = 5;  % MV
        o.m0.sur.pointData.class(ids>cnp(5) & ids<=cnp(6)) = 6;  % TV         
        o.m0.sur.pointData.class(ids>cnp(6) & ids<=cnp(7)) = 7;  % Av
        o.m0.sur.pointData.class(ids>cnp(7)) = 8;  % Pv
    else 
        o.cfg.CobivecoX = false;
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
    end
    mappingDist = sqrt(sum((o.m0.sur.points-surAppended.points(ids,:)).^2,2));
    o.m0.sur.pointData.class(mappingDist > o.cfg.mappingTol*o.m0.meanEdgLen) = 0;
end

P0 = double(o.m0.vol.points);
C0 = double(o.m0.vol.cells);
o.m0.L = cotmatrix(P0, C0);

% If Cobiveco is run, truncSeptSur needs to have only two elements
if o.cfg.CobivecoX == false
    o.cfg.truncSeptSur = o.cfg.truncSeptSur(1:end-1);
end

if o.cfg.exportLevel > 1
    o.m0.debug = o.m0.vol;
    o.m0.debug.pointData.surClass = repmat(uint8(0),size(o.m0.debug.points,1),1);
    o.m0.debug.pointData.surClass(o.m0.surToVol) = o.m0.sur.pointData.class;
end

o.printStatus(sprintf('%.1f seconds\n', toc-t), true);
o.available.mesh0 = true;

end