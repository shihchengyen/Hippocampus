% Run this function with a list of days to classify all cells in the given
% sessions. Outputs a txt file containing the classification of the given
% cell, as well as the p-values of all significance tests run for different
% models on the cell.

function hardcastle_days(tbin_size, fc, day_list)
    % PARAMETERS:
    % tbin_size - size of time bin (in seconds) for binning of vmpv data.
    % fc - number of folds for cross-validation in glm_hardcastle.
    % day_list - array of days to process (as either numbers or strings)
    
    curr_dir = pwd;
    
    for i = 1:length(day_list)
        
        day_dir = day_list(i);
        if ~isa(day_dir, 'char')
            day_dir = num2str(day_dir);
        end
        
        cd(day_dir);
        
        hardcastle_session(tbin_size, fc);
        
    end
    
    cd(curr_dir);

end
