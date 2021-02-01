addpath('..');

%% Create cobiveco object and compute coordinates

c = cobiveco(struct('inPrefix','geo2/heart', 'outPrefix','result_geo2/'));
c.computeAll;
