% Run this function with a list of days, while in directory containing all
% day folders, to classify all cells in the given days. Outputs a txt file 
% containing the classification of the given cell, as well as the p-values 
% of all significance tests run for different models on the cell.

function hardcastle_days(tbin_size, fc, day_list, redo)
    % PARAMETERS:
    % tbin_size - size of time bin (in seconds) for binning of vmpv data.
    % fc - number of folds for cross-validation in glm_hardcastle.
    % day_list - array of days to process (as either numbers or strings)
    
    if ~exist('redo', 'var')
        redo = false;
    end
    
    curr_dir = pwd;
    
    for i = 1:length(day_list)
        % Convert day directory to string if it is not already
        day_dir = day_list(i);
        if ~isa(day_dir, 'char')
            day_dir = num2str(day_dir);
        end
        % Enter day directory
        cd([curr_dir, '/', day_dir]);
        % Process all cells in the session
        hardcastle_session(tbin_size, fc, redo);
        
    end
    
    cd(curr_dir);

end
