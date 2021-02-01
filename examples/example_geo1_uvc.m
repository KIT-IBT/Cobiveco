addpath('..');

%%
c = uvc(struct('inPrefix','geo1/heart_for_UVC', 'outPrefix','result_geo1/', 'exportLevel',2));
c.computeAll;
