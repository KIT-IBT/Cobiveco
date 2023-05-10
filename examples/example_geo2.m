addpath('..');

%% Create cobiveco object and compute coordinates

c = cobiveco(struct('inPrefix','geo2/heart', 'outPrefix','result_geo2/'));

c.prepareMesh0;
if c.cfg.CobivecoX == true
    c.computeAllCobivecoX;
else 
    c.computeAllCobiveco;
end