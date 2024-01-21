function computeApicobasalBridges(o)

% Computes the FreeToSeptalWall (fsw) coordinate on the given object,
% i.e. the coordinate describing the distance from the free wall (1)
% to the septal wall (1.5) on the intervalvular regions.
% The coordinate value is added in the field "ab" to the output objects.
%
% computeApicobasalBridges(o)
%
% Input: 
%   o: instance of the class cobiveco
% 
% Outputs: 
%   o.m0RvBridge, object of class cobiveco [struct] (For details see cobiveco class documentation)
%   o.m0LvBridge, object of class cobiveco [struct] (For details see cobiveco class documentation)
%  

    if ~o.available.rotationalbridges
        o.computeRotationalBridges;
    end

    o.printStatus('Creating Apicobasal on bridges.');
    t = toc;

    o.printStatus('Starting with  RV.');
    % calculate apicobasal coordinate on rv bridge
    
    % set tolerance normGradField
    tolGradField = 1e-8;

    % create ab coordinate
    tmGrad = normalizedGradField(o.m0RvBridge.G, o.m0RvBridge.vol.pointData.tm, tolGradField, true, o.m0RvBridge.vol.points, o.m0RvBridge.vol.cells);
    rtLaplaceGrad = normalizedGradField(o.m0RvBridge.G, o.rtLaplaceBridgeRv, tolGradField, true,  o.m0RvBridge.vol.points, o.m0RvBridge.vol.cells);

    abGrad = cross(tmGrad, rtLaplaceGrad);

    % define boundaries for apicobasal
    idsBoundaryAb1Rv = o.m0RvBridge.surToVol(o.m0RvBridge.sur.pointData.class == 3);
    idsBoundaryAb2Rv = o.m0RvBridge.surToVol(o.m0RvBridge.sur.pointData.class == 4);

    dLvBridgePostAb = solveTrajectDist(o.GrvBridge, abGrad, idsBoundaryAb1Rv, zeros(size(idsBoundaryAb1Rv)), o.cfg.tol, o.cfg.maxit);
    dLvBridgeAntAb = solveTrajectDist(o.GrvBridge, -abGrad, idsBoundaryAb2Rv, zeros(size(idsBoundaryAb2Rv)), o.cfg.tol, o.cfg.maxit);
    abTrajectDistSept = dLvBridgePostAb./(dLvBridgeAntAb+dLvBridgePostAb);

    % Scale coordinate fsw
    abSept = 1 + 0.5 * abTrajectDistSept;
    %fix any outliers
    fixInds = abSept>1.5 | abSept < 0;
    if any(fixInds)
        warning('Needing to fix some fsw for the RV bridge')
        %set them as the mean of their neighbors
        perNeighborCells = o.m0RvBridge.vol.cells(fixInds,:);
        PerNeighborValues = abSept(perNeighborCells);
        neighborsToExclude = intersect(find(fixInds),unique(perNeighborCells(:)));
        PerNeighborValues(ismember(perNeighborCells,neighborsToExclude)) = nan;
        abSept(fixInds) = mean(PerNeighborValues,2,'omitnan');
        %Barycentric interp from good nodes to bad
        %need to first clip out bad nodes
%         perElementKeepMask = ~fixInds(o.m0RvBridge.vol.cells);
%         keepMask_elements = sum(perElementKeepMask,2) >= 4;
%         [keepNodeInds,~,reindexVals] = unique(o.m0RvBridge.vol.cells(keepMask_elements,:)');
%         otherNodeInds = setdiff([1:length(fixInds)],keepNodeInds);
%         newNodeInds = 1:length(keepNodeInds);
%         keepElements = o.m0RvBridge.vol.cells(keepMask_elements,:);
%         clippedElements = reshape(newNodeInds(reindexVals),size(keepElements));
%         M = baryInterpMat(o.m0RvBridge.vol.points(keepNodeInds,:), clippedElements, o.m0RvBridge.vol.points(otherNodeInds,:));
%         abSept(otherNodeInds) = M* abSept(keepNodeInds);
    end

    o.m0RvBridge.vol.cellData.abGrad = single(abGrad);
    o.m0RvBridge.vol.pointData.drvbridgePost = single(dLvBridgePostAb);
    o.m0RvBridge.vol.pointData.drvbridgeAnt = single(dLvBridgeAntAb);
    o.m0RvBridge.vol.pointData.ab = single(abSept);
    vtkWrite(o.m0RvBridge.vol, sprintf('%srvBridgeRtAb.vtu', o.cfg.outPrefix));

    o.printStatus('Continuing with LV.');
    tmGrad = normalizedGradField(o.m0LvBridge.G, o.m0LvBridge.vol.pointData.tm, tolGradField, true, o.m0LvBridge.vol.points, o.m0LvBridge.vol.cells);
    rtLaplaceGrad = normalizedGradField(o.m0LvBridge.G, o.rtLaplaceBridgeLv, tolGradField, true,  o.m0LvBridge.vol.points, o.m0LvBridge.vol.cells);

    abGrad = cross(tmGrad, rtLaplaceGrad);

    % set boundaries for ab coordinate
    idsBoundaryAb1Lv = o.m0LvBridge.surToVol(o.m0LvBridge.sur.pointData.class == 3);
    idsBoundaryAb2Lv = o.m0LvBridge.surToVol(o.m0LvBridge.sur.pointData.class == 4);

    dLvBridgePostAb = solveTrajectDist(o.GlvBridge, abGrad, idsBoundaryAb1Lv, zeros(size(idsBoundaryAb1Lv)), o.cfg.tol, o.cfg.maxit);
    dLvBridgeAntAb = solveTrajectDist(o.GlvBridge, -abGrad, idsBoundaryAb2Lv, zeros(size(idsBoundaryAb2Lv)), o.cfg.tol, o.cfg.maxit);
    abTrajectDistSept = dLvBridgePostAb./(dLvBridgeAntAb+dLvBridgePostAb);

    % Scale coordinate fsw
    abSept = 1 + 0.5 * abTrajectDistSept;
    %fix any outliers
    fixInds = abSept>1.5 | abSept < 0;
    if any(fixInds)
        perNeighborCells = o.m0LvBridge.vol.cells(fixInds,:);
        PerNeighborValues = abSept(perNeighborCells);
        neighborsToExclude = intersect(find(fixInds),unique(perNeighborCells(:)));
        PerNeighborValues(ismember(perNeighborCells,neighborsToExclude)) = nan;
        abSept(fixInds) = mean(PerNeighborValues,2,'omitnan');
    end

    o.m0LvBridge.vol.cellData.abGrad = single(abGrad);
    o.m0LvBridge.vol.pointData.dlvbridgePost = single(dLvBridgePostAb);
    o.m0LvBridge.vol.pointData.dlvbridgeAnt = single(dLvBridgeAntAb);
    o.m0LvBridge.vol.pointData.ab = single(abSept);
    vtkWrite(o.m0LvBridge.vol, sprintf('%slvBridgeRtAb.vtu', o.cfg.outPrefix));

    o.printStatus(sprintf('%.1f seconds\n', toc-t), true);
    o.available.apicobasalbridges = true;


    end