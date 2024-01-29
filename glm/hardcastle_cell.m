% Run this function from cell directory to classify the given cell.
% Outputs a txt file containing the classification of the given
% cell, as well as the p-values of all significance tests run for different
% models on the cell.

function hardcastle_cell(tbin_size, fc, redo)
    % PARAMETERS:
    % tbin_size - size of time bin (in seconds) for binning of vmpv data.
    % fc - number of folds for cross-validation in glm_hardcastle.
    % redo - true/1 or false/0, whether to redo model fitting even if data
    % already exists
    
    if ~exist('redo', 'var')
        redo = false;
    end

    save_dir = [num2str(tbin_size*1000), 'ms_', num2str(fc), 'fold'];

    if ~redo && exist([save_dir, '/glm_hardcastle_results.mat'], 'file')
       % Load pre-generated data file
       load([save_dir, '/glm_hardcastle_results.mat'], 'hc_results');
       if (hc_results.tbin_size ~= tbin_size) || (hc_results.num_folds ~= fc)
           % Need to rerun with corret parameters
           redo = true;
       end
    else
       redo = true;
    end

    if redo
       if exist([save_dir, '/vmpvData.mat'], 'file')
           % Load pre-generated data file
           load([save_dir, '/vmpvData.mat'], 'vmpvData');
           if vmpvData.tbin_size ~= tbin_size
               % Regenerate vmpvData file with correct tbin_size
               vmpvData = glm_vmpvData(tbin_size);
           end
       else
           % Generate vmpvData file from vmpv object
           vmpvData = glm_vmpvData(tbin_size);
       end
       % Run LNP model fitting
       hc_results = glm_hardcastle(vmpvData, fc);
    end

    % Run model forward selection with significance level of 0.05
    [selected_model, hc_results] = select_best_model(hc_results, 0.05);

    % Compare fitted model to ratemaps from actual data using cosine
    % similarity score
    [~, hc_results] = compare_ratemaps(hc_results, 'cosine_similarity');

    % Move results into a labelled directory
    if ~exist(save_dir, 'dir')
       mkdir(save_dir);
    end
    if exist('glm_hardcastle_results.mat', 'file')
       movefile('glm_hardcastle_results.mat', save_dir);
    end
    if exist('vmpvData.mat', 'file')
       movefile('vmpvData.mat', save_dir);
    end
    cd(save_dir);

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

end
