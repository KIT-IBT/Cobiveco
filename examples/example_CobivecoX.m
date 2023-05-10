addpath('..');

%% Create cobiveco object providing all parameters

% config.inPrefix = 'cap1/heart';
% config.outPrefix = 'result_cap1/';
% config.exportLevel = 3;
% config.verbose = true;
% config.mappingTol = 0.5;
% config.tol = 1e-8;
% config.maxit = 3000;
% config.truncSeptSur = [20 10];
% config.abExtrapSmooth = 0.25;
% config.mmgSizingParam = [0.1 0.9 1.1];
% c = cobiveco(config);

%% Create cobiveco object providing only subset of parameters

% c = cobiveco(struct('inPrefix','geo2/geo2', 'outPrefix','result_geo2/', 'exportLevel',3));
c = cobiveco(struct('inPrefix','meanTOF/meanTOF', 'outPrefix','result_meanTOF/', 'exportLevel',3, 'maxit', 4000));

%% Compute coordinates and retrieve config and results

% prepareMesh0 decides on which Cobiveco version to be used
c.prepareMesh0;
if c.cfg.CobivecoX == true
    c.computeAllCobivecoX;
else 
    c.computeAllCobiveco;
end
config = c.cfg;
result = c.result;
R = c.R;

