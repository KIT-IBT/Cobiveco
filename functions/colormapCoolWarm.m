function map = colormapCoolWarm(m)
    % Produces a colormap with m colors that is equivalent to the "Cool to Warm" map in Paraview
    % See https://www.kennethmoreland.com/color-maps/, DOI: 10.1007/978-3-642-10520-3_9
    
    s = linspace(0,1,m);
    rgb1 = [0.230, 0.299, 0.754];
    rgb2 = [0.706, 0.016, 0.150];
    map = colormapDiverging(s, rgb1, rgb2);
end