% Universal Ventricular Coordinates (according to Bayer 2018)
%
% Usage:
%  c = uvc(config);
%  c.computeAll;
%  result = c.result;
%
% Input:
%  - config:                configuration struct with the following fields:
%    config.inPrefix:       prefix for input files, e.g. 'inputDir/namePrefix';
%                           a tetrahedral mesh has to be provided as VTU (inPrefix + '.vtu');
%                           boundaries have to be provided as a VTP file (inPrefix + '.vtp') 
%                           containing point data 'class' with the following values:
%                           1: base, 2: epi, 3: endo_lv, 4: endo_rv, 5: apex, 6: septum_rv
%    config.outPrefix:      prefix for output files, e.g. 'outputDir/' (trailing slash for directories!)
%    config.exportLevel:    integer defining what should be exported as files:
%                           0: nothing
%                           1: only the final result and the rotation matrix (default)
%                           2: + intermediate PDE-solutions, all volume meshes and mmg command line outputs
%                           3: + septal edges, heart axes, apicobasal streamlines and an axes-aligned version of the final result
%    config.verbose:        whether to print status messages (optional, default: true)
%    config.mappingTol:     tolerance for mapping of boundary surfaces, unit: mean edge length (optional, default: 0.5)
%    config.tol:            tolerance for iterative solvers (optional, default: 1e-8)
%    config.maxit:          maximum number of iterations for iterative solvers (optional, default: 3000)
%
% Outputs:
%  - c.result:                 VTK struct with the following point data fields:
%    c.result.pointData.tv:    transventricular coordinate
%    c.result.pointData.tm:    transmural coordinate
%    c.result.pointData.ab:    apicobasal coordinate
%    c.result.pointData.rt:    rotational coordinate
%    c.result.pointData.rtSin: rotational sine coordinate
%    c.result.pointData.rtCos: rotational cosine coordinate
%  - c.cfg:                    complete configuration struct with all parameters
%
% Dependencies:
%  - gptoolbox by Alec Jacobson
%    https://github.com/alecjacobson/gptoolbox
%    Tested version: Commit bee2edb, 20 Aug 2020
%  - vtkToolbox by Steffen Schuler
%    https://github.com/KIT-IBT/vtkToolbox
%  Run dependencies/install.sh to install them from source

classdef uvc < handle
    
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
            'maxit', 3000 ...
            );
        
        % both ventricles
        bv = struct( ...
            'vol', [], ...
            'sur', [], ...
            'surToVol', [], ...
            'L', [], ...
            'tv', [], ...
            'tm', [], ...
            'rtSin', [], ...
            'rtCos', [], ...
            'rt', [], ...
            'abLaplace', [], ...
            'ab', [] ...
            );
            
        % left ventricle
        lv = struct( ...
            'vol', [], ...
            'ids', [], ...
            'sur', [], ...
            'surToVol', [], ...
            'L', [], ...
            'tm', [], ...
            'lap1', [], ...
            'lap2', [], ...
            'lap3', [], ...
            'lap4', [], ...
            'rt', [] ...
            );
            
        % right ventricle
        rv = struct( ...
            'vol', [], ...
            'ids', [], ...
            'sur', [], ...
            'surToVol', [], ...
            'L', [], ...
            'tm', [], ...
            'rt', [] ...
            );
        
        % VTK struct with final coordinates as pointData
        result
        
    end
    
    properties (Access = private)
       
        meanEdgLen
        lvApexVec
        rvApexVec
        
        % used to keep track of what has already been computed
        available = struct( ...
            'mesh', false, ...
            'transventricular', false, ...
            'splitMesh', false, ...
            'transmural', false, ...
            'apicobasal', false, ...
            'rotational', false ...
            );
        
    end

    methods
        
        function o = uvc(config)
            % c = uvc(config)
            
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
            if ~exist('perform_fast_marching_mesh.m', 'file')
                if ~exist([mpath '/../dependencies/numerical-tours'], 'dir')
                    error('Dependency ''numerical-tours'' not found. Run dependencies/install_uvc.sh to install.');
                end
                addpath([mpath '/../dependencies/numerical-tours/matlab/toolbox_general']);
                addpath([mpath '/../dependencies/numerical-tours/matlab/toolbox_graph']);
            end
            if ~exist('cotmatrix.m', 'file')
                if ~exist([mpath '/../dependencies/gptoolbox'], 'dir')
                    error('Dependency ''gptoolbox'' not found. Run dependencies/install_uvc.sh to install.');
                end
                addpath([mpath '/../dependencies/gptoolbox/matrix']);
                addpath([mpath '/../dependencies/gptoolbox/mesh']);
            end
            if ~exist('vtkRead.m', 'file')
                if ~exist([mpath '/../dependencies/vtkToolbox'], 'dir')
                    error('Dependency ''vtkToolbox'' not found. Run dependencies/install_uvc.sh to install.');
                end
                addpath([mpath '/../dependencies/vtkToolbox/MATLAB']);
            end
        end
        
        function computeAll(o)
            % Compute all coordinates
            
            t = toc;
            o.prepareMesh;
            o.computeTransventricular;
            o.splitMesh;
            o.computeTransmural;
            o.computeApicobasal;
            o.computeRotational;
            o.printStatus('Total elapsed time:');
            o.printStatus(sprintf('%.1f seconds\n', toc-t), true);
        end
        
        % declare methods defined in separate files
        prepareMesh(o)
        computeTransventricular(o)
        splitMesh(o)
        computeTransmural(o)
        computeApicobasal(o)
        computeRotational(o)
        
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