function computeApicobasalVentricles(o)

% Computes the apicobasal coordinate on the given object (ventricles),
% i.e. the coordinate refering to the distance from apex (0) to base (1).
% The coordinate value is added in the field "ab" to the output objects.
%
% computeApicobasalVentricles(o)
%
% Input:
%   o: instance of the class cobiveco
%
% Outputs:
%   o.m1, object of class cobiveco [struct] (For details see cobiveco class documentation)
%   o.m0, object of class cobiveco [struct] (For details see cobiveco class documentation)
%

if ~o.available.rotationalbridges
    o.computeRotationalBridges;
end

o.printStatus('Computing Apicobasal coordinate on ventricles.');
t = toc;

% parameters
numTmVal = 10; % number of contour lines in transmural direction
numRtVal = 48; % number of contour lines in rotational direction, ideally a multiple of 6
numClPoints = 100; % number of points per contour line after resampling using cubic splines

% compute contour surfaces of tm (cs)
% for now o.m1Ventricles equals o.m1.vol
tmp = vtkDeleteDataArrays(o.m1Ventricles.vol);
tmp.pointData.tm = o.m1Ventricles.tm;
tmVal = linspace(0.5/numTmVal, 1-0.5/numTmVal, numTmVal);
cs = vtkContourFilter(tmp, 'points', 'tm', tmVal);

% assign regions to cs
cs = vtkConnectivityFilter(cs);
cs.pointData.csRegion = cs.pointData.RegionId;
cs = vtkDeleteDataArrays(cs, {'pointData.csRegion'});

% interpolate abLaplace from o.m2.vol ( o.m2Ventricles) to cs
Mcs = baryInterpMat(o.m2Ventricles.vol.points, o.m2Ventricles.vol.cells, cs.points);
cs.pointData.abLaplace = Mcs * o.m2Ventricles.abLaplace;

if o.cfg.exportLevel > 2
    vtkWrite(cs, sprintf('%sabContourSurfaces.vtp', o.cfg.outPrefix));
end

% find apex point for each csRegion
csRegions = unique(cs.pointData.csRegion);
apexPoints = NaN(numel(csRegions),3);
for i = 1:numel(csRegions)
    pointIds = find(cs.pointData.csRegion == csRegions(i));
    [~,minId] = min(cs.pointData.abLaplace(pointIds));
    apexPoints(i,:) = cs.points(pointIds(minId),:);
end

% compute contour lines of rtSin (rtSinCl)
tmp = cs;
tmp.pointData.rtSin = Mcs * o.m2Ventricles.rtSin;
rtSinVal = sin(pi*linspace(-0.25, 0.25, round(numRtVal/2)+1));
rtSinCl = vtkContourFilter(tmp, 'points', 'rtSin', rtSinVal);
rtSinCl.pointData = rmfield(rtSinCl.pointData, 'rtSin');

% compute contour lines of rtCos (rtCosCl)
tmp = cs;
tmp.pointData.rtCos = Mcs * o.m2Ventricles.rtCos;
rtCosVal = rtSinVal(2:end-1);
rtCosCl = vtkContourFilter(tmp, 'points', 'rtCos', rtCosVal);
rtCosCl.pointData = rmfield(rtCosCl.pointData, 'rtCos');

% combine rtSin and rtCos contour lines (cl)
cl = vtkAppendPolyData({rtSinCl, rtCosCl});

% exclude points close to the apex to split the contour lines into two parts
% apexRadius needs to be large enough for each contour line to be splitted
apexRadius = 3*o.m1Ventricles.meanEdgLen;
cl.pointData.clRegion = zeros(size(cl.points,1),1);
for i = 1:numel(csRegions)
    pointIds = find(cl.pointData.csRegion == csRegions(i));
    ids = rangesearch(cl.points(pointIds,:), apexPoints(i,:), apexRadius);
    cl.pointData.clRegion(pointIds(ids{1})) = -1;
end
cl = vtkThreshold(cl, 'points', 'clRegion', [0 inf]);

% identify separate contour lines --> clRegion
cl = vtkConnectivityFilter(cl);
cl.pointData.clRegion = cl.pointData.RegionId;
cl.pointData = rmfield(cl.pointData, 'RegionId');
cl = rmfield(cl, 'cellData');

% exclude too short contour lines based on values of abLaplace
clRegions = unique(cl.pointData.clRegion);
abLaplaceMin = NaN(numel(clRegions),1);
abLaplaceMax = NaN(numel(clRegions),1);
pointIds = cell(numel(clRegions),1);
for i = 1:numel(clRegions)
    pointIds{i} = find(cl.pointData.clRegion == clRegions(i));
    abLaplaceVal = cl.pointData.abLaplace(pointIds{i});
    abLaplaceMin(i) = min(abLaplaceVal);
    abLaplaceMax(i) = max(abLaplaceVal);
end
abLaplaceMinThresh = min(abLaplaceMin) + 2*(median(abLaplaceMin)-min(abLaplaceMin));
abLaplaceMaxThresh = 0.99;
for i = 1:numel(clRegions)
    if abLaplaceMin(i) > abLaplaceMinThresh || abLaplaceMax(i) < abLaplaceMaxThresh
        cl.pointData.clRegion(pointIds{i}) = -1;
    end
end
cl = vtkThreshold(cl, 'points', 'clRegion', [0 inf]);

% reassign clRegion
cl = vtkConnectivityFilter(cl);
cl.pointData.clRegion = cl.pointData.RegionId;
cl.pointData = rmfield(cl.pointData, 'RegionId');
cl = rmfield(cl, 'cellData');


vtkWrite(cl, sprintf('%sabContourLines.vtp', o.cfg.outPrefix));
% use cubic smoothing spline to smooth and resample the contour lines
% and to compute the normalized distance along the contour lines (abSmoothCl)
clRegions = unique(cl.pointData.clRegion);
splineMisfit = 5e-4 * norm(mean(double(o.m1Ventricles.sur.points(o.m1Ventricles.sur.pointData.class == 1,:)),1)-mean(apexPoints,1));
smoothClPoints = NaN(numel(clRegions)*numClPoints, 3);
abSmoothCl = zeros(numel(clRegions)*numClPoints, 1);
abLength = abSmoothCl;
for i = 1:numel(clRegions)
    pointIds = find(cl.pointData.clRegion == clRegions(i));
    [~,sortInd] = sort(cl.pointData.abLaplace(pointIds));
    apexPoint = apexPoints(csRegions == cl.pointData.csRegion(pointIds(sortInd(1))),:);
    P = [apexPoint; double(cl.points(pointIds(sortInd),:))];
    d = sqrt(sum(diff(P,1,1).^2,2));
    ind = d < 1e-4*o.m1Ventricles.meanEdgLen;
    P(ind,:) = [];
    d(ind) = [];
    d = [0; cumsum(d)];
    w = [100; ones(size(P,1)-1,1)];
    tol = numel(d)*splineMisfit^2;
    sp = spaps(d, P', tol, w);
    P = fnval(sp, linspace(min(d), max(d), numClPoints))';
    d = [0; cumsum(sqrt(sum(diff(P,1,1).^2,2)))];
    ind = (i-1)*numClPoints+1:i*numClPoints;
    smoothClPoints(ind,:) = P;
    abSmoothCl(ind) = d/max(d);
    abLength(ind) = max(d);
end

clSmooth.points = single(smoothClPoints);
clSmooth.pointData.ab = single(abSmoothCl);
clSmooth.cells = int32(repmat(reshape(repmat(0:numClPoints:numel(clRegions)*numClPoints-1,numClPoints-1,1),[],1),1,2) + repmat([(1:numClPoints-1)' (2:numClPoints)'],numel(clRegions),1));
clSmooth.cellTypes = repmat(uint8(3), size(clSmooth.cells,1), 1);
vtkWrite(clSmooth, sprintf('%sabContourLinesSmooth.vtp', o.cfg.outPrefix));

% BEGIN: Laplacian extrapolation of abSmoothCl to o.m1.vol

% ||M ab - abSmoothCl|| forces extrapolated values to fit to contour values

M = baryInterpMat(o.m1Ventricles.vol.points, o.m1Ventricles.vol.cells, smoothClPoints);

% ||E ab - 1|| forces extrapolated values to 1 at the base
% mapping using NN
o.m1Ventricles.sur = vtkArrayMapperNearestNeighbor(o.m2Ventricles.sur, o.m1Ventricles.sur);
o.m1Ventricles.vol.pointData.newsurClass = repmat(uint8(0),size(o.m1Ventricles.vol.points,1),1);
o.m1Ventricles.vol.pointData.newsurClass(o.m1Ventricles.surToVol) = o.m1Ventricles.sur.pointData.newclass;
vtkWrite(o.m1Ventricles.vol, sprintf('%snewmesventricularmeshRemeshedClassestransfered.vtu', o.cfg.outPrefix));
% add new surfaces
baseIds = double(o.m1Ventricles.surToVol(o.m1Ventricles.sur.pointData.newclass == 5 |o.m1Ventricles.sur.pointData.newclass == 6 | o.m1Ventricles.sur.pointData.newclass == 7| o.m1Ventricles.sur.pointData.newclass == 8| o.m1Ventricles.sur.pointData.newclass == 9| o.m1Ventricles.sur.pointData.newclass == 10));
baseVal = ones(numel(baseIds),1);
E = sparse(1:numel(baseIds), baseIds, ones(size(baseIds)), numel(baseIds), size(o.m1Ventricles.vol.points,1));
eta = (numel(abSmoothCl)/(numel(baseIds)))^2;

% ||L ab|| forces extrapolated values to be smooth
L = o.m1Ventricles.massMat \ o.m1Ventricles.L;

% Initial guess for ab using nearest-neighbor interpolation
ab = abSmoothCl(knnsearch(smoothClPoints, o.m1Ventricles.vol.points, 'NSMethod','kdtree'));

% Initial guess for lambda based on the relation between lambda and the
% half-width at half-maximum of the point spread function of the operator
% inv(speye(size(L))+lambda*(L'*L))
hwhm = mean(abLength)/100;
lambda = 1.58*hwhm^3.57;
lambda = numel(abSmoothCl)/numel(ab)*lambda;
fprintf('\n0\t%.3e\n', lambda);

b = M'*abSmoothCl + eta*E'*baseVal;
MM = M'*M + eta*(E'*E);
LL = L'*L;
extrapMisfit = o.cfg.abExtrapSmooth/100;

try
    [~,flag,iter] = secant(@objFun, 1e-1*lambda, lambda, 1e-2*extrapMisfit);
    if flag
        warning('Secant stopped at iteration %i without converging.', iter);
    end
catch err
    disp(getReport(err, 'extended', 'hyperlinks', 'off'));
    warning('Optimization of lambda failed. Using default value instead.');
    objFun(lambda);
end


function objVal = objFun(lambda)
    A = MM + lambda*LL;
    icMat = ichol_autocomp(A, struct('michol','on'));
    [ab, flag, relres, iter] = pcg(A, b, o.cfg.tol, o.cfg.maxit, icMat, icMat', ab);
    if flag
        error('pcg failed at iteration %i with flag %i and relative residual %.1e.', iter, flag, relres);
    end
    objVal = rms(M*ab-abSmoothCl)-extrapMisfit;
end

% END: Laplacian extrapolation

o.m1Ventricles.ab = min(max(ab,0),1);
o.m1Ventricles.debug.pointData.ab = single(o.m1Ventricles.ab);
if o.cfg.exportLevel > 2
    vtkWrite(o.m1Ventricles.debug, sprintf('%so.m1VentriclesFinalAb.vtu', o.cfg.outPrefix));
end

% interpolate to original mesh
P2 = double(o.m1Ventricles.vol.points);
C2 = double(o.m1Ventricles.vol.cells);
o.m1Ventricles.L = cotmatrix(P2, C2);
o.m1Ventricles.G = grad(P2, C2);

idxLvBridge = find(o.m0.vol.pointData.RegionIdVentricle == 2);
idxRvBridge = find(o.m0.vol.pointData.RegionIdVentricle == 1);
idxVentricle = find(o.m0.vol.pointData.RegionIdVentricle == 0);
% map to ventricle
M = baryInterpMat(P2, C2, o.m0.vol.points);
m0AbVentrcile = min(max(M*o.m1Ventricles.ab,0),1);

% map to rv intervalvular region
M2 = baryInterpMat(o.m0RvBridge.vol.points,o.m0RvBridge.vol.cells, o.m0.vol.points);
abRv = M2*double(o.m0RvBridge.vol.pointData.ab);

% map to lv intervalvular region
M3 = baryInterpMat(o.m0LvBridge.vol.points,o.m0LvBridge.vol.cells, o.m0.vol.points);
abLv = M3*double(o.m0LvBridge.vol.pointData.ab);

o.m0.vol.pointData.ab = repmat(uint8(0),size(o.m0.vol.points,1),1);
o.m0.vol.pointData.ab = m0AbVentrcile;
o.m0.vol.pointData.ab(idxLvBridge) = abLv(idxLvBridge);
o.m0.vol.pointData.ab(idxRvBridge) = abRv(idxRvBridge);
vtkWrite(o.m0.vol, sprintf('%soriginalmeshAbScaled.vtu', o.cfg.outPrefix));

if o.cfg.exportLevel > 2
    o.m0.vol.pointData.bridgefilter = repmat(uint8(0),size(o.m0.vol.points,1),1);
    o.m0.vol.pointData.bridgefilter(idxVentricle) = 0;
    o.m0.vol.pointData.bridgefilter(idxLvBridge) = 3;
    o.m0.vol.pointData.bridgefilter(idxRvBridge) = 2;
    vtkWrite(o.m0.vol, sprintf('%sbridgeclasses.vtu', o.cfg.outPrefix));

end

% interpolate rotational
% map to ventricle
Mrt = baryInterpMat(o.m2Ventricles.vol.points, o.m2Ventricles.vol.cells, o.m0.vol.points);
m0RtVentricle = min(max(Mrt*o.m2Ventricles.rt,0),1);

% map to rv intervalvular region
M2rt = baryInterpMat(o.m0RvBridge.vol.points,o.m0RvBridge.vol.cells, o.m0.vol.points);
rtRv = M2rt*double(o.m0RvBridge.vol.pointData.rt);

% map to lv intervalvular region
M3rt = baryInterpMat(o.m0LvBridge.vol.points,o.m0LvBridge.vol.cells, o.m0.vol.points);
rtLv = M3rt*double(o.m0LvBridge.vol.pointData.rt);

o.m0.vol.pointData.rt = repmat(uint8(0),size(o.m0.vol.points,1),1);
o.m0.vol.pointData.rt = m0RtVentricle;
o.m0.vol.pointData.rt(idxLvBridge) = rtLv(idxLvBridge);
o.m0.vol.pointData.rt(idxRvBridge) = rtRv(idxRvBridge);

o.m0.vol.pointData.rtSin = sin(2*pi*o.m0.vol.pointData.rt);
o.m0.vol.pointData.rtCos = cos(2*pi*o.m0.vol.pointData.rt);
vtkWrite(o.m0.vol, sprintf('%soriginalmeshRtScaled.vtu', o.cfg.outPrefix));


% add to original mesh
o.m0.rtSin = min(max(o.m2Ventricles.M*o.m2Ventricles.rtSin,-1),1);
o.m0.rtCos = min(max(o.m2Ventricles.M*o.m2Ventricles.rtCos,-1),1);
o.m0.rt = atan2(o.m0.rtSin,o.m0.rtCos)/(2*pi);
o.m0.rt(o.m0.rt<0) = o.m0.rt(o.m0.rt<0)+1;
o.result.pointData.rtSin = o.m0.rtSin;
o.result.pointData.rtCos = o.m0.rtCos;
o.result.pointData.rt = o.m0.rt;
vtkWrite(o.result, sprintf('%sresultRt.vtu', o.cfg.outPrefix));

% as we do not have a discontinuity within intervalvular regions, we can
% map without using the sin and cosine
o.m0LvBridge.rtSin = min(max(M3rt*o.m0LvBridge.vol.pointData.rtSin,-1),1);
o.m0LvBridge.rtCos = min(max(M3rt*o.m0LvBridge.vol.pointData.rtCos,-1),1);
o.m0RvBridge.rtSin = min(max(M2rt*o.m0RvBridge.vol.pointData.rtSin,-1),1);
o.m0RvBridge.rtCos = min(max(M2rt*o.m0RvBridge.vol.pointData.rtCos,-1),1);

o.result.pointData.rtSin = single(o.m0.rtSin);
o.result.pointData.rtSin(idxLvBridge) = o.m0LvBridge.rtSin(idxLvBridge);
o.result.pointData.rtSin(idxRvBridge) = o.m0RvBridge.rtSin(idxRvBridge);
o.result.pointData.rtCos = single(o.m0.rtCos);
o.result.pointData.rtCos(idxLvBridge) = o.m0LvBridge.rtCos(idxLvBridge);
o.result.pointData.rtCos(idxRvBridge) = o.m0RvBridge.rtCos(idxRvBridge);
o.result.pointData.rt(idxRvBridge) = rtRv(idxRvBridge);
o.result.pointData.rt(idxLvBridge) = rtLv(idxLvBridge);

o.result.pointData.ab = single(o.m0.vol.pointData.ab);
o.result.pointData.tv = single(o.m0.tv);
o.result.pointData.tm = single(o.m0.tm);
vtkWrite(o.result, sprintf('%sresultRtWithBridges.vtu', o.cfg.outPrefix));

% rtSin and rtCos are not required for mapping in bridge regions but are 
% currently calculated
o.m0.vol.pointData.rtSin = min(max(Mrt*o.m2Ventricles.rtSin,-1),1);
o.m0.vol.pointData.rtCos = min(max(Mrt*o.m2Ventricles.rtCos,-1),1);
% map from rv
rtRvSin = min(max(M2rt*o.m0RvBridge.vol.pointData.rtSin,-1),1);
rtRvCos = min(max(M2rt*o.m0RvBridge.vol.pointData.rtCos,-1),1);
% map from lv
rtLvSin = min(max(M3rt*o.m0LvBridge.vol.pointData.rtSin,-1),1);
rtLvCos = min(max(M3rt*o.m0LvBridge.vol.pointData.rtCos,-1),1);

% update rsin and rcos with bridge data
o.m0.vol.pointData.rtSin(idxRvBridge) = rtRvSin(idxRvBridge);
o.m0.vol.pointData.rtSin(idxLvBridge) = rtLvSin(idxLvBridge);
o.m0.vol.pointData.rtCos(idxRvBridge) = rtRvCos(idxRvBridge);
o.m0.vol.pointData.rtCos(idxLvBridge) = rtLvCos(idxLvBridge);

% calculate rt
o.m0.vol.pointData.rt = atan2(o.m0.vol.pointData.rtSin,o.m0.vol.pointData.rtCos)/(2*pi);
o.m0.vol.pointData.rt(o.m0.vol.pointData.rt<0) = o.m0.vol.pointData.rt(o.m0.vol.pointData.rt<0)+1;

o.m0.vol.pointData.rtSin = single(o.m0.vol.pointData.rtSin);
o.m0.vol.pointData.rtCos = single(o.m0.vol.pointData.rtCos);
o.m0.vol.pointData.rt = single(o.m0.vol.pointData.rt);


vtkWrite(o.m0.vol, sprintf('%soriginalmeshRt.vtu', o.cfg.outPrefix));

o.printStatus(sprintf('%.1f seconds\n', toc-t), true);
o.available.apicobasalventricles = true;

end
