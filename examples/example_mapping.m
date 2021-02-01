addpath('../utilities');
addpath('../functions');
addpath('HeatGeodesic');

%% Define paths

sourcePrefix = 'result_geo1/';
targetPrefix = 'result_geo2/';
outputPrefix = 'result_mapping/';

outPath = fileparts(outputPrefix);
if ~isempty(outPath) && ~exist(outPath, 'dir')
    mkdir(outPath);
end

%% Load geometries

source = vtkRead(sprintf('%sresult.vtu', sourcePrefix));
target = vtkRead(sprintf('%sresult.vtu', targetPrefix));

%% Compute mapping matrix

searchradius = 2;
M = cobiveco_computeMappingMatrix(source, target, 'linear', searchradius, true);

%% Align geometries with global axes
%  This is optional for mapping data (example 1 below),
%  but required for obtaining a mean geometry (example 2 below)

R_source = load(sprintf('%sR.mat', sourcePrefix)); R_source = R_source.R;
R_target = load(sprintf('%sR.mat', targetPrefix)); R_target = R_target.R;
source = cobiveco_applyAlignmentMatrix(source, R_source);
target = cobiveco_applyAlignmentMatrix(target, R_target);

%% Example 1: Mapping of a geodesic distance field

P = double(source.points);
C = double(source.cells);
L = cotmatrix(P,C);
G = grad(P,C);
D = div(P,C);

t = 2;
tol = 1e-8;
maxit = 1000;
id = 15926;

d = computeHeatGeodesic(L, G, D, id, t, tol, maxit);

src = vtkDeleteDataArrays(source);
src.pointData.d = d;
vtkWrite(src, sprintf('%ssrc.vtu', outputPrefix));

tar = vtkDeleteDataArrays(target);
tar.pointData.d = M*d;
vtkWrite(tar, sprintf('%star.vtu', outputPrefix))

%% Example 2: Mapping of node coords to obtain a mean geometry

mean = target;
mean.points = single((mean.points + M*double(source.points))/2);
vtkWrite(mean, sprintf('%smean.vtu', outputPrefix))
