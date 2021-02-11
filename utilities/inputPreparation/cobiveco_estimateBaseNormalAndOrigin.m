function [baseNormal,baseOrigin,debug] = cobiveco_estimateBaseNormalAndOrigin(sur)

for k = 1:5

    if k == 1
        axSur = sur;
        baseNormal = [1 0 0];
    else
        tSur = vtkConnectivityFilter(vtkThreshold(sur, 'points', 'height', [prctile(sur.pointData.height,p) inf]));
        reg = cell(3,1);
        r = NaN(3,1);
        reg{1} = vtkThreshold(tSur, 'points', 'RegionId', [0 0]);
        reg{2} = vtkThreshold(tSur, 'points', 'RegionId', [1 1]);
        reg{3} = vtkThreshold(tSur, 'points', 'RegionId', [2 2]);
        [~,r(1)] = vtkSmallestEnclosingSphere(reg{1});
        [~,r(2)] = vtkSmallestEnclosingSphere(reg{2});
        [~,r(3)] = vtkSmallestEnclosingSphere(reg{3});
        [~,epiInd] = max(r);
        endoInd = setdiff(1:3,epiInd);
        axSur = vtkAppendPolyData({reg{endoInd(1)}, reg{endoInd(2)}});
    end
    baseNormal = computeLongAxis(axSur, baseNormal);
    sur.pointData.height = double(sur.points*baseNormal');
    
    p = 50; % initial percentile value used for thresholding
    for i = 1:5
        tSur = vtkConnectivityFilter(vtkThreshold(sur, 'points', 'height', [prctile(sur.pointData.height,p) inf]));
        if numel(unique(tSur.pointData.RegionId)) < 3
            if k == 1 && i == 1
                tSur = vtkConnectivityFilter(vtkThreshold(sur, 'points', 'height', [-inf prctile(sur.pointData.height,p)]));
                if numel(unique(tSur.pointData.RegionId)) < 3
                    error('There must be 3 regions after initial thresholding.');
                end
                baseNormal = -baseNormal;
                sur.pointData.height = -sur.pointData.height;
            else
                break;
            end
        end
        p = p-10;
    end
    p = p+9;
    for i = 1:9
        tSur = vtkConnectivityFilter(vtkThreshold(sur, 'points', 'height', [prctile(sur.pointData.height,p) inf]));
        if numel(unique(tSur.pointData.RegionId)) < 3
            break;
        end
        p = p-1;
    end
    p = p+1;

end

baseHeight = prctile(sur.pointData.height,p);
center = vtkCenterOfArea(axSur);
centerHeight = center*baseNormal';
baseOrigin = center + (baseHeight-centerHeight)*baseNormal;
sur.pointData.height = sur.pointData.height-baseHeight;

centerDist = double((sur.points-center)*baseNormal');
axis = vtkCreateStruct([center+1.2*min(centerDist)*baseNormal; center; center+1.2*max(centerDist)*baseNormal], [1 2; 2 3]);
axis.pointData.height = double(axis.points*baseNormal'-baseHeight);

debug = vtkAppendPolyData({sur, axis});

end