% This example calculates the two way mapping error as described by Bayer et al., 2018 by 
% mapping CobivecoX coordinates between source and target.
% Here, the user can add an excel file (ID) containing the names of the 
% targets to map coordinates to and from.

addpath('../utilities');
addpath('../functions');
addpath('../dependencies/vtkToolbox')


% set up search paths
mpath = fileparts(mfilename('fullpath'));
addpath([mpath '/../functions']);
addpath('../utilities');
addpath('../utilities/inputPreparation');
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

%% Define data paths
ID = readmatrix('List_of_models_to_map.xlsx', 'OutputType','string');
sourcePrefix = './result_source/';
originalPrefix = './original/';
source = vtkRead(sprintf('%sresult.vtu', sourcePrefix));

original_euclidean_coords = vtkRead(sprintf('%ssource.vtu', originalPrefix));
euclidean_original = original_euclidean_coords.points;
source.pointData.euclidean_coords = euclidean_original;

% split into bridges and ventricles
s_ab = source.pointData.ab;
s_ventricles = find(s_ab<=1);
s_bridges = find(s_ab>1);
source.pointData.ventricles = zeros(size(source.pointData.ab,1),1);
source.pointData.ventricles(s_ventricles) = 1;
bridge_struct_src = source;
ventricle_src = source;
ventricle_src_celldata = vtkPointDataToCellData(ventricle_src);
bridge_struct_src_celldata = vtkPointDataToCellData(bridge_struct_src);
s_ventricles = vtkDeleteDataArrays(vtkThreshold(ventricle_src_celldata, 'cells', 'ventricles', [1 1]));
s_bridges = vtkDeleteDataArrays(vtkThreshold(bridge_struct_src_celldata, 'cells', 'ventricles', [0 0]));

ventricle_errors = zeros(size(s_ventricles.points,1),1);
bridge_errors = zeros(size(s_bridges.points,1),1);

for i=1:size(ID,1)
    try
        fprintf('Working on ID #%s\n', ID(i));
        targetPrefix = sprintf('./result_%s/', ID(i));
        outputPrefix = sprintf('./result_mapping_two_way_%s/', ID(i));
        outPath = fileparts(outputPrefix);

        if ~isempty(outPath) && ~exist(outPath, 'dir')
            mkdir(outPath);
        end

        target = vtkRead(sprintf('%sresult.vtu', targetPrefix));

        % separate into ventricles and bridges
        t_ab = target.pointData.ab;
        s_ab = source.pointData.ab;
        s_ventricles = find(s_ab<=1);
        s_bridges = find(s_ab>1);
        t_ventricles = find(t_ab<=1);
        t_bridges = find(t_ab>1);
        source.pointData.ventricles = zeros(size(source.pointData.ab,1),1);
        source.pointData.ventricles(s_ventricles) = 1;
        target.pointData.ventricles = zeros(size(target.pointData.ab,1),1);
        target.pointData.ventricles(t_ventricles) = 1;
        bridge_struct_src = source;
        ventricle_src = source;
        bridge_struct_tar = target;
        ventricle_tar = target;
        ventricle_src_celldata = vtkPointDataToCellData(ventricle_src);
        target_celldata = vtkPointDataToCellData(target);
        bridge_struct_src_celldata = vtkPointDataToCellData(bridge_struct_src);
        ventricle_tar_celldata = vtkPointDataToCellData(ventricle_tar);

        s_ventricles = vtkDeleteDataArrays(vtkThreshold(ventricle_src_celldata, 'cells', 'ventricles', [1 1]));
        t_ventricles = vtkDeleteDataArrays(vtkThreshold(target_celldata, 'cells', 'ventricles', [1 1]));
        s_bridges = vtkDeleteDataArrays(vtkThreshold(bridge_struct_src_celldata, 'cells', 'ventricles', [0 0]));
        t_bridges = vtkDeleteDataArrays(vtkThreshold(ventricle_tar_celldata, 'cells', 'ventricles', [0 0]));

        % map data
        s_ventricles = vtkArrayMapperNearestNeighbor(source, s_ventricles);
        t_ventricles = vtkArrayMapperNearestNeighbor(target, t_ventricles);
        t_bridges = vtkArrayMapperNearestNeighbor(target, t_bridges);
        s_bridges = vtkArrayMapperNearestNeighbor(source, s_bridges);
        vtkWrite(t_ventricles, sprintf('%star_ventricles.vtu', outputPrefix))
        vtkWrite(s_ventricles, sprintf('%ss_ventricles.vtu', outputPrefix))
        vtkWrite(t_bridges, sprintf('%st_bridges.vtu', outputPrefix))
        vtkWrite(s_bridges, sprintf('%ss_bridges.vtu', outputPrefix))
        
        %% Compute mapping matrix
        searchradius = 3;
        M = cobiveco_computeMappingMatrix(s_ventricles, t_ventricles, 'linear', searchradius, true);
        t_ventricles.pointData.euclidean_coords = M*s_ventricles.pointData.euclidean_coords;
        %% map back
        searchradius = 3;
        M_back = cobiveco_computeMappingMatrix(t_ventricles, s_ventricles, 'linear', searchradius, true);
        s_ventricles.pointData.mapped_euclidean = M_back*t_ventricles.pointData.euclidean_coords;
        EuclideanDistance = vecnorm(s_ventricles.pointData.mapped_euclidean-s_ventricles.pointData.euclidean_coords, 2,2);
        s_ventricles.pointData.EuclideanDistance = EuclideanDistance;

        %% Map Bridges

        %% Compute mapping matrix
        searchradius = 3;
        M = cobiveco_computeMappingMatrix_briges(s_bridges, t_bridges, 'linear', searchradius, true);
        t_bridges.pointData.euclidean_coords = M*s_bridges.pointData.euclidean_coords;

        %% map back
        % Compute mapping matrix
        searchradius = 3;
        M_back = cobiveco_computeMappingMatrix_briges(t_bridges, s_bridges, 'linear', searchradius, true);
        s_bridges.pointData.mapped_euclidean = M_back*t_bridges.pointData.euclidean_coords;
        EuclideanDistance_b = vecnorm(s_bridges.pointData.mapped_euclidean-s_bridges.pointData.euclidean_coords, 2,2);
        s_bridges.pointData.EuclideanDistance = EuclideanDistance_b;

        % add new result to arrays
        if i == 1
            ventricle_errors = EuclideanDistance;
            bridge_errors = EuclideanDistance_b;
        else
            ventricle_errors = [ventricle_errors EuclideanDistance];
            bridge_errors = [bridge_errors EuclideanDistance_b];
        end

        fprintf('Finished ID #%s\n', ID(i));
        fid = fopen('Finished_IDs','a+');
        fprintf(fid, 'instance %s \n', ID(i));
        fclose(fid);
    catch ME
        disp('Error Message:')
        disp(ME.message)
        fid = fopen('Failed_IDs','a+');
        fprintf(fid, 'Failed instance %s\n', ID(i));
        fprintf(fid, '%s', ME.getReport('extended', 'hyperlinks','off'));
        fclose(fid);
        fprintf('Failed on ID #%s\n', ID(i));
        continue;
    end

end


total_ventricle_error = ventricle_errors;
total_bridge_error = bridge_errors;
% calculate the mean over all mappings
mean_ventricles = mean(total_ventricle_error,2);
mean_bridge = mean(total_bridge_error,2);
% get the total mean
mean_bridge = mean(total_bridge_error,"all");
mean_ventricles = mean(total_ventricle_error,'all');

std_ventricles = std(total_ventricle_error, 1, 'all');
std_bridges = std(total_bridge_error, 1, 'all');

number_mappings = size(total_ventricle_error,2);

fid = fopen('Mapping output.txt', 'wt');
fprintf(fid, 'mean_bridge = %d \n', mean_bridge);
fprintf(fid, 'mean_ventricles = %d \n', mean_ventricles);
fprintf(fid, 'std_bridges = %d \n', std_bridges);
fprintf(fid, 'std_ventricles = %d \n', std_ventricles);
fprintf(fid, 'number mappings = %d \n', number_mappings);
fprintf(fid, 'ids mapped  %d \n', number_mappings);
fprintf(fid, 'instance%s \n', ID);
fclose(fid);