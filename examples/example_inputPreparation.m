addpath('../utilities');
addpath('../utilities/inputPreparation');
addpath('../functions');

%% Read original mesh

vol = vtkRead('geo3/heart_original.vtu');

%% Estimate the normal vector and the origin of a basal plane

sur = vtkDataSetSurfaceFilter(vol);
[baseNormal,baseOrigin,debug] = cobiveco_estimateBaseNormalAndOrigin(sur);
vtkWrite(debug, 'geo3/debug1.vtp');

%% Adjust baseNormal and baseOrigin, if needed

baseShift = -7;
%baseShift = -15;
baseOrigin = baseOrigin + baseShift*baseNormal;

%% Clip mesh at the basal plane

vol = cobiveco_clipBase(vol, baseNormal, baseOrigin);
vtkWrite(vol, 'geo3/debug2.vtu');

%% Create surface classes

sur = vtkDataSetSurfaceFilter(vol);
maxAngle = 40; % max angle of face normals wrt baseNormal for defining the base class
numSubdiv = 2; % can help for coarse meshes (interpolation of face normals)
[sur,debug] = cobiveco_createClasses(sur, baseNormal, maxAngle, numSubdiv);
vtkWrite(debug, 'geo3/debug3.vtp');

%% Remove bridges
%keyboard
[vol,debug,mmgOutput] = cobiveco_removeBridges(vol, sur, baseNormal, 'rv', true);
vtkWrite(debug, 'geo3/debug4.vtp');

%% Recreate surface classes
%keyboard
sur = vtkDataSetSurfaceFilter(vol);
[sur,debug] = cobiveco_createClasses(sur, baseNormal, maxAngle, numSubdiv);
vtkWrite(debug, 'geo3/debug5.vtp');

%% Write result

vtkWrite(sur, 'geo3/test_clipping_heart.vtp');
vtkWrite(vol, 'geo3/test_clipping_heart.vtu');
