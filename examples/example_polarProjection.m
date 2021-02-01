addpath('../utilities');
addpath('../functions');
addpath('HeatGeodesic');

%% Load geometry and create some data to plot

source = vtkRead('result_geo1/result.vtu');
data = source.points(:,1);

%% Precompute polar projection

tm = -1; % transmural coordinate (0: epi, 1: endo, -1: epi and endo)
[M,mask] = cobiveco_createPolarProjection(source, tm, [], [], [], [], [], [], [], true);

%% Create polar projection plot using precomputed M and mask

cobiveco_createPolarProjection(M, mask, data);
