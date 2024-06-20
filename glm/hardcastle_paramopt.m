% Run this function from directory containing vmpv.mat.
% Outputs a txt file containing the classification of the simulated
% cell, as well as the p-values of all significance tests run for different
% models on the cell.

function hardcastle_paramopt(tbin_size, fc, param_range, redo)
    % PARAMETERS:
    % tbin_size - size of time bin (in seconds) for binning of vmpv data.
    % fc - number of folds for cross-validation in glm_hardcastle.
    % param_range - range of beta values to use for testing. Will use the
    % same beta value for all three variables (place, headdirection, spatialview).
    % redo - true/1 or false/0, whether to redo model fitting even if data
    % already exists
   
    if ~exist('param_range', 'var')
        param_range = 3;  % use default smoothing value of 3
    end
   
    if ~exist('redo', 'var')
        redo = false;
    end

    if exist('genData.mat', 'file')
       % Load pre-generated data file
       load('genData.mat', 'genData');
    else
       % Generate simulated data from vmpv object
       genData = glm_genData_gaussian(tbin_size);
    end

    for i = 1:length(param_range)
       param_val = param_range(i);
       save_dir = ['smoothing_', num2str(param_val)];
       smooth_params = [ param_val, param_val, param_val ];

       disp(['Running classification for smoothing beta value of ', num2str(param_val)]);
       mkdir(save_dir);
       cd(save_dir);

       % Run LNP model fitting and model forward selection
       if exist('glm_hardcastle_results.mat', 'file') && ~redo
          load('glm_hardcastle_results.mat', 'hc_results');
       else
           hc_results = glm_hardcastle(genData, fc, smooth_params);
       end
       [selected_model, hc_results] = select_best_model(hc_results);

       % Cell classification and significance testing, saved as a txt file
       % in the labelled directory
       if exist('hardcastle_class.txt', 'file')
           % Overwrite old text file
           delete('hardcastle_class.txt')
       end
       diary('hardcastle_class.txt');
       disp(['Cell classification: ' selected_model]);
       disp(' ');
       hardcastle_testing_LLH(hc_results);
       diary off;

       % Generate response plots for classified variables
       glm_hardcastle_plot(hc_results, selected_model, [], true);

       cd ..;
    end

    % Save parameters used in the testing run as a txt file
    if exist('run_params.txt', 'file')
       % Overwrite old text file
       delete('run_params.txt')
    end
    diary('run_params.txt');
    disp(['tbin_size: ', num2str(tbin_size)]);
    disp(['fc: ', num2str(fc)]);
    disp(['cell_type: ', genData.cell_params.cell_type]);
    diary off;

end
