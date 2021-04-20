% Consistent Biventricular Coordinates
%
% Usage:
%  c = cobiveco(config);
%  c.computeAll;
%  result = c.result;
%
% Input:
%  - config:                configuration struct with the following fields:
%    config.inPrefix:       prefix for input files, e.g. 'inputDir/namePrefix';
%                           a tetrahedral mesh has to be provided as VTU (inPrefix + '.vtu');
%                           boundary surfaces can either be providing as individual PLY files:
%                           inPrefix + '_base.ply', '_epi.ply', '_endo_lv.ply', '_endo_rv.ply'
%                           or as one VTP file (inPrefix + '.vtp') containing point data 'class' with the following values:
%                           1: base, 2: epi, 3: endo_lv, 4: endo_rv
%    config.outPrefix:      prefix for output files, e.g. 'outputDir/' (trailing slash for directories!)
%    config.exportLevel:    integer defining what should be exported as files:
%                           0: nothing
%                           1: only the final result and the rotation matrix (default)
%                           2: + intermediate PDE-solutions, all volume meshes and mmg command line outputs
%                           3: + septal curves, heart axes, apicobasal streamlines and an axes-aligned version of the final result
%    config.verbose:        whether to print status messages (optional, default: true)
%    config.mappingTol:     tolerance for mapping of boundary surfaces, unit: mean edge length (optional, default: 0.5)
%    config.tol:            tolerance for iterative solvers (optional, default: 1e-8)
%    config.maxit:          maximum number of iterations for iterative solvers (optional, default: 3000)
%    config.truncSeptSur:   percentages of triangles to truncate from the septal surface on the anterior and posterior side, respectively;
%                           affects the left-right axis and the position of the apex point (optional, default: [20 10])
%    config.abExtrapSmooth: smoothness for extrapolation of apicobasal coord from contour lines to mesh nodes
%                           (RMS deviation from the contour line values in percent), useful range: 0.1...1 (optional, default: 0.25)
%    config.mmgSizingParam: a vector [hausd hmin hmax] defining sizing parameters for mmg (optional), unit: mean edge length
%                           hausd: maximal Hausdorff distance for approximation of boundaries (default: 0.1)
%                           hmin:  minimal edge size (default: 0.9)
%                           hmax:  maximal edge size (default: 1.1)
%
% Outputs:
%  - c.result:                 VTK struct with the following point data fields:
%    c.result.pointData.tv:    transventricular coordinate
%    c.result.pointData.tm:    transmural coordinate
%    c.result.pointData.ab:    apicobasal coordinate
%    c.result.pointData.rt:    rotational coordinate
%    c.result.pointData.rtSin: rotational sine coordinate
%    c.result.pointData.rtCos: rotational cosine coordinate
%  - c.R:                      rotation matrix for alignment of heart axes with x,y,z
%  - c.cfg:                    complete configuration struct with all parameters
%
% Dependencies:
%  - MMG3D by INP/Inria/U Bordeaux
%    https://github.com/MmgTools/mmg
%    Tested version: Release 5.4.3, 26 Feb 2020
%  - gptoolbox by Alec Jacobson
%    https://github.com/alecjacobson/gptoolbox
%    Tested version: Commit bee2edb, 20 Aug 2020
%  - vtkToolbox by Steffen Schuler
%    https://github.com/KIT-IBT/vtkToolbox
%  Run dependencies/install.sh to install them from source
%
% Written in 2020 by Steffen Schuler
% Institute of Biomedical Engineering, KIT
% www.ibt.kit.edu

classdef cobiveco < handle
    
    properties (SetAccess = private)
        
        % Configuration parameters
        % (default values)
        cfg = struct( ...
            'inPrefix', '', ...
            'outPrefix', '', ...
            'exportLevel', 1, ...
            'verbose', true, ...
            'mappingTol', 0.5, ...
            'tol', 1e-8, ...
            'maxit', 3000, ...
            'truncSeptSur', [20 10], ...
            'abExtrapSmooth', 0.25, ...
            'mmgSizingParam', [0.1 0.9 1.1] ...
            );
        
        % Mesh 0 (original mesh)
        m0 = struct( ...
            'vol', [], ...
            'sur', [], ...
            'surToVol', [], ...
            'meanEdgLen', [], ...
            'L', [], ...
            'tv', [], ...
            'tm', [], ...
            'rtSin', [], ...
            'rtCos', [], ...
            'rt', [], ...
            'ab', [], ...
            'debug', [] ...
            );
        
        % Mesh 1 (isovalue discretization at boundary between LV and RV)
        m1 = struct( ...
            'vol', [], ...
            'sur', [], ...
            'surToVol', [], ...
            'meanEdgLen', [], ...
            'L', [], ...
            'G', [], ...
            'massMat', [], ...
            'M', [], ...
            'tm', [], ...
            'ridgeLaplace', [], ...
            'debug', [] ...
            );
        
        % Mesh 2 (isovalue discretization at ridge)
        m2 = struct( ...
            'vol', [], ...
            'sur', [], ...
            'surToVol', [], ...
            'L', [], ...
            'G', [], ...
            'M', [], ...
            'abLaplace', [], ...
            'rtSin', [], ...
            'rtCos', [], ...
            'rt', [], ...
            'debug', [] ...
            );
        
        % VTK struct with final coordinates as pointData
        result
        
        % Rotation matrix aligning heart axes with x,y,z
        R
        
    end
    
    properties (Access = private)
       
        septCurveAnt
        septCurvePost
        
        % used to keep track of what has already been computed
        available = struct( ...
            'mesh0', false, ...
            'transventricular', false, ...
            'mesh1', false, ...
            'transmural', false, ...
            'heartAxesAndApex', false, ...
            'mesh2', false, ...
            'rotational', false, ...
            'apicobasal', false ...
            );
        
    end

    methods
        
        function o = cobiveco(config)
            % c = cobiveco(config)
            
            % overwrite default parameters, if provided by the user
            fns = fieldnames(config);
            for i = 1:numel(fns)
                fn = fns{i};
                if ~isfield(o.cfg,fn)
                    error('Unknown parameter ''%s'' found in config.', fn);
                end
                o.cfg.(fn) = config.(fn);
            end
            
            % check mandatory parameters
            if ~isfield(config,'inPrefix')
                error('Mandatory parameter ''inPrefix'' missing in config.');
            end
            if o.cfg.exportLevel > 0
                if isempty(o.cfg.outPrefix)
                    error('Mandatory parameter ''outPrefix'' missing in config or empty.');
                end
                outPath = fileparts(o.cfg.outPrefix);
                if ~isempty(outPath) && ~exist(outPath, 'dir')
                    mkdir(outPath);
                end
                testfile = sprintf('%stest', o.cfg.outPrefix);
                [fid,errmsg] = fopen(testfile, 'w');
                if fid==-1
                    error('Invalid parameter outPrefix. Could not create ''%s'' (%s).', testfile, errmsg);
                else
                    fclose(fid);
                    delete(testfile);
                end
            end
            
            % start timer
            tic;
            
            % set up paths
            mpath = fileparts(mfilename('fullpath'));
            addpath([mpath '/../functions']);
            if ~exist([mpath '/../dependencies/mmg/build/bin/mmg3d_O3'], 'file')
                error('Dependency ''mmg'' not found. Run dependencies/install_cobiveco.sh to install.');
            end
            if ~exist('cotmatrix.m', 'file')
                if ~exist([mpath '/../dependencies/gptoolbox'], 'dir')
                    error('Dependency ''gptoolbox'' not found. Run dependencies/install_cobiveco.sh to install.');
                end
                addpath([mpath '/../dependencies/gptoolbox/matrix']);
                addpath([mpath '/../dependencies/gptoolbox/mesh']);
            end
            if ~exist('vtkRead.m', 'file')
                if ~exist([mpath '/../dependencies/vtkToolbox'], 'dir')
                    error('Dependency ''vtkToolbox'' not found. Run dependencies/install_cobiveco.sh to install.');
                end
                addpath([mpath '/../dependencies/vtkToolbox/MATLAB']);
            end
        end
        
        function computeAll(o)
            % Compute all coordinates
            
            t = toc;
            o.prepareMesh0;
            o.computeTransventricular;
            o.prepareMesh1;
            o.computeTransmural;
            o.computeHeartAxesAndApex;
            o.prepareMesh2;
            o.computeRotational;
            o.computeApicobasal;
            o.exportResult;
            o.printStatus('Total elapsed time:');
            o.printStatus(sprintf('%.1f seconds\n', toc-t), true);
        end
        
        % declare methods defined in separate files
        prepareMesh0(o)
        computeTransventricular(o)
        prepareMesh1(o)
        computeTransmural(o)
        computeHeartAxesAndApex(o)
        prepareMesh2(o)
        computeRotational(o)
        computeApicobasal(o)
        exportResult(o)
        
    end
    
    methods (Access = private)
        
        function printStatus(o, msg, noPadding)
            if o.cfg.verbose
                if nargin > 2 && noPadding
                    fprintf('%s', msg);
                else
                    pad = repmat(' ', 1, 41-numel(msg));
                    fprintf('%s%s', msg, pad);
                end
            end
        end
        
    end
    
end