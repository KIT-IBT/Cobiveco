
addpath('..');

%% Create cobiveco object providing all parameters

% config.inPrefix = 'directory/heart';
% config.outPrefix = 'result_heart/';
% config.exportLevel = 3;
% config.verbose = true;
% config.mappingTol = 0.5;
% config.tol = 1e-8;
% config.maxit = 3000;
% config.truncSeptSur = [20 10];
% config.abExtrapSmooth = 0.25;
% config.mmgSizingParam = [0.1 0.9 1.1];
% c = cobiveco(config);

% read in list of strings representing names of your models and directories
ID = readmatrix('Some_excel_file.xlsx', OutputType='string');
% or add array
%ID = ["some_name";"some_other_name"; "yet_another_name"];

for i=1:size(ID,1)
    try
        %% Create cobiveco object providing only subset of parameters
        c = cobiveco(struct('inPrefix',sprintf('%s/%s',ID(i),ID(i)), 'outPrefix',sprintf('result_%s/',ID(i)), 'exportLevel',3, 'tol', 1e-8, 'maxit', 4000));
        
        %% Compute coordinates and retrieve config and results
        c.prepareMesh0;
        if c.cfg.CobivecoX == true
            c.computeAllCobivecoX;
        else 
            c.computeAllCobiveco;
        end
        config = c.cfg;
        result = c.result;
        R = c.R;
        fprintf('Finished ID #%s\n', ID(i));
    catch ME
        disp('Error Message:')
        disp(ME.message)
        fid = fopen('List_of_errors','a+');
        fprintf(fid, '---------------------------------- \n');
        fprintf(fid, 'Failed instance %s \n', ID(i));
        fprintf(fid, '%s \n', ME.getReport('extended', 'hyperlinks','off'));
        fprintf(fid, '---------------------------------- \n');
            fclose(fid);
        fprintf('Failed on ID #%s\n', ID(i));
    continue;
    end

end

