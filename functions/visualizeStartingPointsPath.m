function visualizeStartingPointsPath(struct1, point1, point2, point3, point4, point5, point6, point7, point8, cfg)

    % Visualize the starting points used for the shortest graph.
    %
    % visualizeStartingPointsPath(struct1, point1, point2, point3, point4, point5, point6, point7, point8, cfg)
    % 
    % Inputs:
    %
    %   struct1 [numPoints x 3]
    %   point1 double, [1x3]
    %   point2 double, [1x3]
    %   point3 double, [1x3]
    %   point4 double, [1x3]
    %   point5 double, [1x3]
    %   point6 double, [1x3]
    %   point7 double, [1x3]
    %   point8 double, [1x3]
    %   cfg, object/structure in cobiveco class which contains info about configuration storage location
    %
    %
    % Written by Lisa Pankewitz
    
    struct1.vol.pointData.startingpoints = zeros(size(struct1.vol.points,1),1);
    point1InGraph = ismembertol(struct1.vol.points,point1,0.1,'ByRows',true,'DataScale',1.0);
    [idxInGraphPoint1,values] = find(point1InGraph==1);
    struct1.vol.pointData.startingpoints(idxInGraphPoint1) = 1;

    point2InGraph = ismembertol(struct1.vol.points,point2,0.1,'ByRows',true,'DataScale',1.0);
    [idxInGraphPoint2,values] = find(point2InGraph==1);
    struct1.vol.pointData.startingpoints(idxInGraphPoint2) = 2;

    point3InGraph = ismembertol(struct1.vol.points,point3,0.1,'ByRows',true,'DataScale',1.0);
    [idxInGraphPoint3,values] = find(point3InGraph==1);
    struct1.vol.pointData.startingpoints(idxInGraphPoint3) = 3;

    point4InGraph = ismembertol(struct1.vol.points,point4,0.1,'ByRows',true,'DataScale',1.0);
    [idxInGraphPoint4,values] = find(point4InGraph==1);
    struct1.vol.pointData.startingpoints(idxInGraphPoint4) = 4;

    point5InGraph = ismembertol(struct1.vol.points,point5,0.1,'ByRows',true,'DataScale',1.0);
    [idxInGraphPoint5,values] = find(point5InGraph==1);
    struct1.vol.pointData.startingpoints(idxInGraphPoint5) = 5;

    point6InGraph = ismembertol(struct1.vol.points,point6,0.1,'ByRows',true,'DataScale',1.0);
    [idxInGraphPoint6,values] = find(point6InGraph==1);
    struct1.vol.pointData.startingpoints(idxInGraphPoint6) = 6;

    point7InGraph = ismembertol(struct1.vol.points,point7,0.1,'ByRows',true,'DataScale',1.0);
    [idxInGraphPoint7,values] = find(point7InGraph==1);
    struct1.vol.pointData.startingpoints(idxInGraphPoint7) = 7;

    point8InGraph = ismembertol(struct1.vol.points,point8,0.1,'ByRows',true,'DataScale',1.0);
    [idxInGraphPoint8,values] = find(point8InGraph==1);
    struct1.vol.pointData.startingpoints(idxInGraphPoint8) = 8;

    vtkWrite(struct1.vol, sprintf('%sstarts.vtu', cfg.outPrefix));

end