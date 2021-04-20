function computeApicobasal(o)

if ~o.available.rotational
    o.computeRotational;
end
o.printStatus('Computing apicobasal coordinate...');
t = toc;

% parameters
numTmVal = 10; % number of contour lines in transmural direction
numRtVal = 48; % number of contour lines in rotational direction, ideally a multiple of 6
numClPoints = 100; % number of points per contour line after resampling using cubic splines

% compute contour surfaces of tm (cs)
tmp = vtkDeleteDataArrays(o.m1.vol);
tmp.pointData.tm = o.m1.tm;
tmVal = linspace(0.5/numTmVal, 1-0.5/numTmVal, numTmVal);
cs = vtkContourFilter(tmp, 'points', 'tm', tmVal);

% assign regions to cs
cs = vtkConnectivityFilter(cs);
cs.pointData.csRegion = cs.pointData.RegionId;
cs = vtkDeleteDataArrays(cs, {'pointData.csRegion'});

% interpolate abLaplace from o.m2.vol to cs
Mcs = baryInterpMat(o.m2.vol.points, o.m2.vol.cells, cs.points);
cs.pointData.abLaplace = Mcs * o.m2.abLaplace;

if o.cfg.exportLevel > 2
    vtkWrite(cs, sprintf('%sabContourSurfaces.vtp', o.cfg.outPrefix));
end

% find apex point for each csRegion
csRegions = unique(cs.pointData.csRegion);
apexPoints = NaN(numel(csRegions),3);
for i = 1:numel(csRegions)
    pointIds = find(cs.pointData.csRegion==csRegions(i));
    [~,minId] = min(cs.pointData.abLaplace(pointIds));
    apexPoints(i,:) = cs.points(pointIds(minId),:);
end

% compute contour lines of rtSin (rtSinCl)
tmp = cs;
tmp.pointData.rtSin = Mcs * o.m2.rtSin;
rtSinVal = sin(pi*linspace(-0.25, 0.25, round(numRtVal/2)+1));
rtSinCl = vtkContourFilter(tmp, 'points', 'rtSin', rtSinVal);
rtSinCl.pointData = rmfield(rtSinCl.pointData, 'rtSin');

% compute contour lines of rtCos (rtCosCl)
tmp = cs;
tmp.pointData.rtCos = Mcs * o.m2.rtCos;
rtCosVal = rtSinVal(2:end-1);
rtCosCl = vtkContourFilter(tmp, 'points', 'rtCos', rtCosVal);
rtCosCl.pointData = rmfield(rtCosCl.pointData, 'rtCos');

% combine rtSin and rtCos contour lines (cl)
cl = vtkAppendPolyData({rtSinCl, rtCosCl});

% exclude points close to the apex to split the contour lines into two parts
% apexRadius needs to be large enough for each contour line to be splitted
apexRadius = 3*o.m1.meanEdgLen;
cl.pointData.clRegion = zeros(size(cl.points,1),1);
for i = 1:numel(csRegions)
    pointIds = find(cl.pointData.csRegion==csRegions(i));
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
    pointIds{i} = find(cl.pointData.clRegion==clRegions(i));
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

if o.cfg.exportLevel > 2
    vtkWrite(cl, sprintf('%sabContourLines.vtp', o.cfg.outPrefix));
end

% use cubic smoothing spline to smooth and resample the contour lines
% and to compute the normalized distance along the contour lines (abSmoothCl)
clRegions = unique(cl.pointData.clRegion);
splineMisfit = 5e-4 * norm(mean(double(o.m1.sur.points(o.m1.sur.pointData.class==1,:)),1)-mean(apexPoints,1));
smoothClPoints = NaN(numel(clRegions)*numClPoints, 3);
abSmoothCl = zeros(numel(clRegions)*numClPoints, 1);
abLength = abSmoothCl;
for i = 1:numel(clRegions)
    pointIds = find(cl.pointData.clRegion==clRegions(i));
    [~,sortInd] = sort(cl.pointData.abLaplace(pointIds));
    apexPoint = apexPoints(csRegions==cl.pointData.csRegion(pointIds(sortInd(1))),:);
    P = [apexPoint; double(cl.points(pointIds(sortInd),:))];
    d = sqrt(sum(diff(P,1,1).^2,2));
    ind = d < 1e-4*o.m1.meanEdgLen;
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

if o.cfg.exportLevel > 2
    clSmooth.points = single(smoothClPoints);
    clSmooth.pointData.ab = single(abSmoothCl);
    clSmooth.cells = int32(repmat(reshape(repmat(0:numClPoints:numel(clRegions)*numClPoints-1,numClPoints-1,1),[],1),1,2) + repmat([(1:numClPoints-1)' (2:numClPoints)'],numel(clRegions),1));
    clSmooth.cellTypes = repmat(uint8(3), size(clSmooth.cells,1), 1);
    vtkWrite(clSmooth, sprintf('%sabContourLinesSmooth.vtp', o.cfg.outPrefix));
end

% BEGIN: Laplacian extrapolation of abSmoothCl to o.m1.vol

% ||M ab - abSmoothCl|| forces extrapolated values to fit to contour values
M = baryInterpMat(o.m1.vol.points, o.m1.vol.cells, smoothClPoints);

% ||E ab - 1|| forces extrapolated values to 1 at the base
baseIds = double(o.m1.surToVol(o.m1.sur.pointData.class==1));
baseVal = ones(numel(baseIds),1);
E = sparse(1:numel(baseIds), baseIds, ones(size(baseIds)), numel(baseIds), size(o.m1.vol.points,1));
eta = (numel(abSmoothCl)/(numel(baseIds)))^2;

% ||L ab|| forces extrapolated values to be smooth
L = o.m1.massMat \ o.m1.L;

% Initial guess for ab using nearest-neighbor interpolation
ab = abSmoothCl(knnsearch(smoothClPoints, o.m1.vol.points, 'NSMethod','kdtree'));

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

o.m1.ab = min(max(ab,0),1);
o.m0.ab = min(max(o.m1.M*o.m1.ab,0),1);
o.result.pointData.ab = single(o.m0.ab);

if o.cfg.exportLevel > 1
    o.m1.debug.pointData.ab = single(o.m1.ab);
    vtkWrite(o.m1.debug, sprintf('%sdebug1.vtu', o.cfg.outPrefix));
    
    o.m0.debug.pointData.ab = single(o.m0.ab);
    vtkWrite(o.m0.debug, sprintf('%sdebug0.vtu', o.cfg.outPrefix));
end

o.printStatus(sprintf('%.1f seconds\n', toc-t), true);
o.available.apicobasal = true;

end
