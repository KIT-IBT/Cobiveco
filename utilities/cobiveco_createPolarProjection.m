function [M,mask,polar] = cobiveco_createPolarProjection(source, tm, data, r, ax, annotateLines, annotateText, method, searchradius, verbose)

if nargin < 10
    verbose = false;
end
if nargin < 9
    searchradius = [];
end
if nargin < 8
    method = [];
end
if nargin < 7 || isempty(annotateText)
    annotateText = true;
end
if nargin < 6 || isempty(annotateLines)
    annotateLines = true;
end
if nargin < 5
    ax = [];
end
if nargin < 4 || isempty(r)
    r = 200;
end

if ~isa(data, 'double')
    data = double(data);
end

%% Parameters for annotations

lineColor = [0.2 0.2 0.2];
markerSize = 20;
fontSize = 16;

%% Set up cartesian grid and transform it into polar coords

if ~issparse(source) && ~islogical(tm)
    precomputed = false;
    [x,y] = meshgrid(-r:r);
    [theta, rho] = cart2pol(x,y);
    mask = rho <= r;
    
elseif issparse(source) && islogical(tm)
    precomputed = true;
    M = source;
    mask = tm;
    switch size(M,1)/nnz(mask)
        case 1 % epicardium
            tm = 0;
        case 2 % endocardium/below-epicardium
            tm = 1;
        case 3 % epi- and endocardium
            tm = -1;
        otherwise
            error('M and mask do not fit to each other.');
    end
    r = (length(mask)-1)/2;
    
else
    error('Both M and mask have been provided instead of source and tm to run in precomputed mode.');
end

%% Compute polar plot(s)

if tm == 0 % epicardium
    
    if ~precomputed
        left = theta>=-pi/2 & theta<=pi/2;
        thetaEpi = theta;
        thetaEpi(left) = 4/3*(thetaEpi(left)+pi/2);
        thetaEpi(1:2*r+1,1:r+1) = fliplr(thetaEpi(1:2*r+1,r+1:2*r+1));

        target.pointData.ab = rho(mask)/r;
        target.pointData.rtSin = sin(thetaEpi(mask));
        target.pointData.rtCos = cos(thetaEpi(mask));
        target.pointData.tm = zeros(nnz(mask),1);
        target.pointData.tv = [ones(floor(nnz(mask)/2),1); zeros(ceil(nnz(mask)/2),1)];

        M = cobiveco_computeMappingMatrix(source, target, method, searchradius, verbose);
    end
    
    if ~isempty(data)
        polar = NaN(2*r+1,2*r+1);
        polar(mask) = M * data;
        polar = polar(2:end-1,2:end-1);
    end
    
elseif tm > 0 && tm <= 1 % endocardium / below-epicardium
    
    if ~precomputed
        target.pointData.ab    = repmat(rho(mask)/r, 2, 1);
        target.pointData.rtSin = [sin(-theta(mask)-1/3*pi); sin(theta(mask)-4/3*pi)];
        target.pointData.rtCos = [cos(-theta(mask)-1/3*pi); cos(theta(mask)-4/3*pi)];
        target.pointData.tm    = repmat(tm, 2*nnz(mask), 1);
        target.pointData.tv    = [ones(nnz(mask),1); zeros(nnz(mask),1)];

        M = cobiveco_computeMappingMatrix(source, target, method, searchradius, verbose);
    end

    if ~isempty(data)
        polar_rv = NaN(2*r+1,2*r+1);
        polar_rv(mask) = M(1:nnz(mask),:) * data;
        polar_rv = polar_rv(2:end-1,2:end-1);

        polar_lv = NaN(2*r+1,2*r+1);
        polar_lv(mask) = M(nnz(mask)+1:end,:) * data;
        polar_lv = polar_lv(2:end-1,2:end-1);

        polar = [polar_rv NaN(2*r-1,1) polar_lv];
    end
    
elseif tm == -1 % epi- and endocardium
    
    if ~precomputed
        left = theta>=-pi/2 & theta<=pi/2;
        thetaEpi = theta;
        thetaEpi(left) = 4/3*(thetaEpi(left)+pi/2);
        thetaEpi(1:2*r+1,1:r+1) = fliplr(thetaEpi(1:2*r+1,r+1:2*r+1));
        
        target.pointData.ab    = repmat(rho(mask)/r, 3, 1);
        target.pointData.rtSin = [sin(thetaEpi(mask)); sin(-theta(mask)-1/3*pi); sin(theta(mask)-4/3*pi)];
        target.pointData.rtCos = [cos(thetaEpi(mask)); cos(-theta(mask)-1/3*pi); cos(theta(mask)-4/3*pi)];
        target.pointData.tm    = [zeros(nnz(mask),1); ones(2*nnz(mask),1)];
        target.pointData.tv    = [ones(floor(nnz(mask)/2),1); zeros(ceil(nnz(mask)/2),1); ones(nnz(mask),1); zeros(nnz(mask),1)];

        M = cobiveco_computeMappingMatrix(source, target, method, searchradius, verbose);
    end

    if ~isempty(data)
        polar_epi = NaN(2*r+1,2*r+1);
        polar_epi(mask) = M(1:nnz(mask),:) * data;
        polar_epi = polar_epi(2:end-1,2:end-1);

        polar_rv = NaN(2*r+1,2*r+1);
        polar_rv(mask) = M(nnz(mask)+1:2*nnz(mask),:) * data;
        polar_rv = polar_rv(2:end-1,2:end-1);

        polar_lv = NaN(2*r+1,2*r+1);
        polar_lv(mask) = M(2*nnz(mask)+1:end,:) * data;
        polar_lv = polar_lv(2:end-1,2:end-1);

        polar = [polar_epi NaN(2*r-1,1) polar_rv NaN(2*r-1,1) polar_lv];
    end
    
else
    error('tm must be in the range [0,1] or -1.')
end

%% Display polar plot(s)

if ~isempty(data)
    
    if isempty(ax) || ~isgraphics(ax)
        figure('DefaultAxesPosition', [0.02, 0.02, 0.96, 0.82]);
        ax = gca;
    end

    imagesc(ax, polar, 'AlphaData',~isnan(polar));
    axis(ax, 'equal');
    xlim(ax, [0 size(M,1)/nnz(mask)*2*r]);
    ylim(ax, [0 2*r]);
    set(ax, 'visible','off');
    colormap(ax, colormapCoolWarm(20));

    %% Add annotations

    hold(ax, 'on');

    if tm == 0 % epicardium

        if annotateLines
            % radial line
            plot(ax, [r r], [0.5 2*r-0.5], 'color',lineColor);

            % circle
            rc = r*cos(0:pi/100:2*pi);
            rs = r*sin(0:pi/100:2*pi);
            plot(ax, r+rc, r+rs, 'color',lineColor);

            % apex point
            plot(ax, r, r, '.', 'MarkerSize',20, 'color',lineColor);
        end

        if annotateText
            % titles
            text(ax, 1/2*r, -r/9, 'RV', 'FontSize',fontSize, 'FontUnits','Normalized', 'HorizontalAlignment','center');
            text(ax, 3/2*r, -r/9, 'LV', 'FontSize',fontSize, 'FontUnits','Normalized', 'HorizontalAlignment','center');
        end

    elseif tm > 0 % endocardium/below-epicardium

        if annotateLines
            % radial lines
            rc = r*cos(pi/3);
            rs = r*sin(pi/3);
            plot(ax, [  r   r+rc], [r r+rs], 'color',lineColor);
            plot(ax, [  r   r+rc], [r r-rs], 'color',lineColor);
            plot(ax, [3*r 3*r-rc], [r r+rs], 'color',lineColor);
            plot(ax, [3*r 3*r-rc], [r r-rs], 'color',lineColor);

            % circles
            rc = r*cos(0:pi/100:2*pi);
            rs = r*sin(0:pi/100:2*pi);
            plot(ax,   r+rc, r+rs, 'color',lineColor);
            plot(ax, 3*r+rc, r+rs, 'color',lineColor);

            % apex points
            plot(ax,   r, r, '.', 'color',lineColor, 'MarkerSize',markerSize);
            plot(ax, 3*r, r, '.', 'color',lineColor, 'MarkerSize',markerSize);
        end

        if annotateText
            % titles
            text(ax,   r, -r/9, 'RV',     'FontSize',fontSize, 'FontUnits','Normalized', 'HorizontalAlignment','center');
            text(ax, 2*r, -r/9, 'Septum', 'FontSize',fontSize, 'FontUnits','Normalized', 'HorizontalAlignment','center');
            text(ax, 3*r, -r/9, 'LV',     'FontSize',fontSize, 'FontUnits','Normalized', 'HorizontalAlignment','center');
        end

    else % epi- and endocardium

        if annotateLines
            % radial lines
            plot(ax, [r r], [0.5 2*r-0.5], 'color',lineColor);
            rc = r*cos(pi/3);
            rs = r*sin(pi/3);
            plot(ax, [3*r 3*r+rc], [r r+rs], 'color',lineColor);
            plot(ax, [3*r 3*r+rc], [r r-rs], 'color',lineColor);
            plot(ax, [5*r 5*r-rc], [r r+rs], 'color',lineColor);
            plot(ax, [5*r 5*r-rc], [r r-rs], 'color',lineColor);

            % circles
            rc = r*cos(0:pi/100:2*pi);
            rs = r*sin(0:pi/100:2*pi);
            plot(ax,   r+rc, r+rs, 'color',lineColor);
            plot(ax, 3*r+rc, r+rs, 'color',lineColor);
            plot(ax, 5*r+rc, r+rs, 'color',lineColor);

            % apex points
            plot(ax,   r, r, '.', 'color',lineColor, 'MarkerSize',markerSize);
            plot(ax, 3*r, r, '.', 'color',lineColor, 'MarkerSize',markerSize);
            plot(ax, 5*r, r, '.', 'color',lineColor, 'MarkerSize',markerSize);
        end

        if annotateText
            % titles
            text(ax,     r, -r/3, 'Epicardium',  'FontSize',fontSize, 'FontWeight','bold', 'FontUnits','Normalized', 'HorizontalAlignment','center');
            text(ax, 1/2*r, -r/9, 'RV',          'FontSize',fontSize, 'FontUnits','Normalized', 'HorizontalAlignment','center');
            text(ax, 3/2*r, -r/9, 'LV',          'FontSize',fontSize, 'FontUnits','Normalized', 'HorizontalAlignment','center');
            text(ax,   4*r, -r/3, 'Endocardium', 'FontSize',fontSize, 'FontWeight','bold', 'FontUnits','Normalized', 'HorizontalAlignment','center');
            text(ax,   3*r, -r/9, 'RV',          'FontSize',fontSize, 'FontUnits','Normalized', 'HorizontalAlignment','center');
            text(ax,   4*r, -r/9, 'Septum',      'FontSize',fontSize, 'FontUnits','Normalized', 'HorizontalAlignment','center');
            text(ax,   5*r, -r/9, 'LV',          'FontSize',fontSize, 'FontUnits','Normalized', 'HorizontalAlignment','center');
        end

    end

    hold(ax, 'off');
    
end

end