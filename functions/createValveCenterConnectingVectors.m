function [TV_PV_vector, PV_TV_vector, MV_AV_vector, AV_MV_vector] = createValveCenterConnectingVectors(center_mv, center_av, center_tv,center_pv)
% Calculate vectors between centers of all valve annuli in biventricular mesh.
%
% [TV_PV_vector, PV_TV_vector, MV_AV_vector, AV_MV_vector] = createValveCenterConnectingVectors(center_mv, ,center_av, center_tv,center_pv)
%
% Inputs:
%   center_mv: center mitrale valve annulus 3d point, double [1x3]
%   center_av: center aortic valve annulus 3d point, double [1x3]
%   center_tv: center tricuspid valve annulus 3d point, double [1x3]
%   center_pv: center pulmonary valve annulus 3d point, double [1x3]
%
% Outputs:
%   Four vectors connecting each of the valves center points: [TV_PV_vector, PV_TV_vector, MV_AV_vector, AV_MV_vector], double [1x3]
%
% Written by Lisa Pankewitz, 2022

% find vector between center TV and PV
TV_PV_vector = center_pv - center_tv;
PV_TV_vector = center_tv - center_pv;

% find vector between center MV and AV
MV_AV_vector = center_av - center_mv;
AV_MV_vector = center_mv - center_av;

end



