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

% read in matrix
%ID = readmatrix('List.xlsx');
ID = [389];

for i=1:size(ID,1)
    % Check if the computer is in sleep mode (only works for Windows)
    [status, result] = system('powercfg /requests');
    if strcmp(result, 'The following sleep states are available: Standby (S1-S3) Hibernate Hybrid Sleep\nThe following sleep states are not available: Standby (S0 Low Power Idle)\n')
      % Pause the Matlab process for 20 minutes
      pause(60 * 20);
    else
        try
            %% Create cobiveco object providing only subset of parameters
            c = cobiveco(struct('inPrefix',sprintf('instance%d/instance%d',ID(i),ID(i)), 'outPrefix',sprintf('result_instance%d/',ID(i)), 'exportLevel',3, 'tol', 1e-8, 'maxit', 4000));

            %% Compute coordinates and retrieve config and results
            fprintf('Working on ID #%d\n', ID(i));
            c.prepareMesh0;
            if c.cfg.CobivecoX == true
                c.computeAllCobivecoX;
            else
                c.computeAllCobiveco;
            end
            config = c.cfg;
            result = c.result;
            R = c.R;
            fprintf('Finished ID #%d\n', ID(i));
        catch ME
            disp('Error Message:')
            disp(ME.message)
            fid = fopen('Error_Log','a+');
            fprintf(fid, '---------------------------------- \n');
            fprintf(fid, 'Failed instance %i \n', ID(i));
            fprintf(fid, '%s \n', ME.getReport('extended', 'hyperlinks','off'));
            fprintf(fid, '---------------------------------- \n');
            fclose(fid);
            fprintf('Failed on ID #%d \n', ID(i));
        continue;% Jump to next iteration of: for i
        end
    end

end

%% Export config as JSON

%fid = fopen(sprintf('%sconfig.json', config.outPrefix), 'w');
%fwrite(fid, jsonencode(config));
%fclose(fid);
