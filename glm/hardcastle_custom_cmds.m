% Custom script to run custom functions on all the cells specified in the
% file given by 'filename'.

function hardcastle_custom_cmds(filename)
    % PARAMETERS:
    % tbin_size - size of time bin (in seconds) for binning of vmpv data.
    % fc - number of folds for cross-validation in glm_hardcastle.
    % cell_list - txt file containing list of cells to process.

    % Read in list of cells from txt file
    if ~exist('filename', 'var')
        filename = 'cell_list.txt';  % should be in the directory level containing all day folders
    end
    cell_list = textread(filename, '%s', 'delimiter', '\n');
    model_types = {{'place', 'headdirection', 'spatialview'}, {'place', 'headdirection'}, {'place', 'spatialview'}, ...
        {'headdirection', 'spatialview'}, {'place'}, {'headdirection'}, {'spatialview'}};
    
    curr_dir = pwd;
    process_dirs = {'50ms_10fold(23-09-18)', '50ms_10fold_nosmooth(23-09-18)'};
    
    for i = 1:length(cell_list)
        % Go into each cell directory, and run hardcastle functions
        cd([curr_dir '/' cell_list{i}]);
        disp(['Processing ', pwd]);
        
        load(['50ms_10fold(23-09-18)' '/glm_hardcastle_results(smoothed).mat'], 'hc_results');
        cell_class = hc_results.classification;
        if isnan(cell_class)
            continue
        end
        var_list = model_types{cell_class};
            
        for k = 1:length(process_dirs)
            dir = process_dirs{k};
            cd(dir);
            
            filenames = {'place_ratemap_fits', 'hd_ratemap_fits', 'view_ratemap_fits'};
            suffixes = {'(raw).png', '(smoothed).png', '(training).png'};
            
            for a = 1:length(filenames)
                for b = 1:length(suffixes)
                    filename = [filenames{a} suffixes{b}];
                    if exist(filename, 'file')
                        if ~exist('full', 'dir')
                            mkdir('full');
                        end
                        movefile(filename, 'full/');
                    end
                end
            end

            for j = 1:length(var_list)
                unix(['/Users/jcheng/miniconda3/bin/python /Users/jcheng/Documents/MATLAB/Hippocampus/glm/plot_ratemap_fits.py', ' ', ...
                    var_list{j}, ' ', 'training', ' ', 'glm_hardcastle_results\(training\).mat']);
                switch var_list{j}
                    case 'place'
                        movefile('place_ratemap_fits.png', 'place_ratemap_fits(training).png')
                    case 'headdirection'
                        movefile('hd_ratemap_fits.png', 'hd_ratemap_fits(training).png')
                    case 'spatialview'
                        movefile('view_ratemap_fits.png', 'view_ratemap_fits(training).png')
                end

                unix(['/Users/jcheng/miniconda3/bin/python /Users/jcheng/Documents/MATLAB/Hippocampus/glm/plot_ratemap_fits.py', ' ', ...
                    var_list{j}, ' ', 'smoothed', ' ', 'glm_hardcastle_results\(smoothed\).mat']);
                switch var_list{j}
                    case 'place'
                        movefile('place_ratemap_fits.png', 'place_ratemap_fits(smoothed).png')
                    case 'headdirection'
                        movefile('hd_ratemap_fits.png', 'hd_ratemap_fits(smoothed).png')
                    case 'spatialview'
                        movefile('view_ratemap_fits.png', 'view_ratemap_fits(smoothed).png')
                end

                unix(['/Users/jcheng/miniconda3/bin/python /Users/jcheng/Documents/MATLAB/Hippocampus/glm/plot_ratemap_fits.py', ' ', ...
                    var_list{j}, ' ', 'raw', ' ', 'glm_hardcastle_results\(raw\).mat']);
                switch var_list{j}
                    case 'place'
                        movefile('place_ratemap_fits.png', 'place_ratemap_fits(raw).png')
                    case 'headdirection'
                        movefile('hd_ratemap_fits.png', 'hd_ratemap_fits(raw).png')
                    case 'spatialview'
                        movefile('view_ratemap_fits.png', 'view_ratemap_fits(raw).png')
                end
            end
            
            cd ..;
        end
    end
    cd(curr_dir);
end
