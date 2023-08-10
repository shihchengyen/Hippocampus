% Run this function from day/session directory to classify all cells in the given
% session. Outputs a txt file containing the classification of the given
% cell, as well as the p-values of all significance tests run for different
% models on the cell.

function hardcastle_batch(tbin_size, fc)
    % PARAMETERS:
    % tbin_size - size of time bin (in seconds) for binning of vmpv data.
    % fc - number of folds for cross-validation in glm_hardcastle.

    % Look for all cells in the current session directory
    [~, cell_list] = unix('find $PWD -name "spiketrain.mat" | rev | cut -d "/" -f 2- | rev');
    cell_list = strsplit(cell_list, '\n');
    cell_list = cell_list(1:end-1); % last entry after strsplit is an empty string so discard it
    
    curr_dir = pwd;
    
    for i = 1:length(cell_list)
        % Go into each cell directory, and run hardcastle functions
       cd(cell_list{i});
       disp(['Running classfication for ', pwd]);

       if exist('glm_hardcastle_results.mat', 'file')
           % Load pre-generated data file
           load('glm_hardcastle_results.mat', 'hc_results');
           if (hc_results.tbin_size ~= tbin_size) || (hc_results.num_folds ~= fc)
               % Regenerate data from vmpv object and run LNP model fitting
               vmpv_data = glm_vmpvData(tbin_size);
               hc_results = glm_hardcastle(vmpv_data, fc);
           end
       else
           % Generate data from vmpv object and run LNP model fitting
           vmpv_data = glm_vmpvData(tbin_size);
           hc_results = glm_hardcastle(vmpv_data, fc);
       end
       % Run model forward selection
       selected_model = select_best_model(hc_results);
       
       % Cell classification and significance testing, saved as a txt file
       % in the cell directory
       if exist('hardcastle_class.txt', 'file')
           % Overwrite old text file
           delete('hardcastle_class.txt')
       end
       diary('hardcastle_class.txt')
       disp(['Cell classification: ' selected_model]);
       disp(' ');
       hardcastle_testing(hc_results);
       diary off
       
       % Generate response plots for classified variables
       glm_hardcastle_plot(hc_results, selected_model);
       
    end
    
    cd(curr_dir);

end
