% Run this function from day/session directory to classify all cells in the given
% session. Outputs a txt file containing the classification of the given
% cell, as well as the p-values of all significance tests run for different
% models on the cell.

function hardcastle_session(tbin_size, fc, redo)
    % PARAMETERS:
    % tbin_size - size of time bin (in seconds) for binning of vmpv data.
    % fc - number of folds for cross-validation in glm_hardcastle.
    
    if ~exist('redo', 'var')
        redo = false;
    end

    % Look for all cells in the current session directory
    [~, cell_list] = unix('find $PWD -name "spiketrain.mat" | rev | cut -d "/" -f 2- | rev');
    cell_list = strsplit(cell_list, '\n');
    cell_list = cell_list(1:end-1); % last entry after strsplit is an empty string so discard it
    
    curr_dir = pwd;
    
    for i = 1:length(cell_list)
        % Go into each cell directory, and run hardcastle functions
       cd(cell_list{i});
       disp(['Running classfication for ', pwd]);
       
       hardcastle_cell(tbin_size, fc, redo);
    
    end
    
    cd(curr_dir);

end
