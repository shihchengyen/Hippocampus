function glm_hardcastle_plot(hc_results, model, axRange, save)
%	Plots fitted params onto each variable space,
%   analogous to the tuning curve of the cell.
%   Plots all folds in a single figure for each model variable.

%	PARAMETERS:
%	hc_results - struct output of glm_hardcastle
%	model - 'place' / 'headdirection' / 'spatialview' / 
%           'ph' / 'pv' / 'hv' / 'phv'
%   save - true/1 or false/0, whether to save the figure(s) as a .fig file

if ~exist('save', 'var')
    save = false;
end
if ~exist('model', 'var')
    model = hc_results.classification;
end
if isnumeric(model)
    model_names = {'phv', 'ph', 'pv', 'hv', 'place', 'headdirection', 'spatialview'};
    if ~isnan(model)
        model = model_names{model};
    else
        model = 'unclassified';
    end
end   

filter_unoccupied = true;  % set unoccupied bins to NaN when plotting
% large_negative_number = -1e3;
large_negative_number = -1e1;
pltRange = [0, 99];  % percentiles to set the colorbar range of the plots to, if not using axRange

params = hc_results.params_consol;
tbin_size = hc_results.tbin_size;
num_folds = hc_results.num_folds;
if isfield(hc_results, 'similarity_scores')
    similarity_scores = hc_results.similarity_scores;
else
    similarity_scores = { nan(num_folds, 1), nan(num_folds, 1), nan(num_folds, 1), ...
        nan(num_folds, 1), nan(num_folds, 1), nan(num_folds, 1), nan(num_folds, 1) };
end

[subplot_rows, subplot_cols] = getSubplotGridSize(num_folds);

%%% Specify environment bin geometry here %%%
% can be changed
floor_width = 40;
wall_height = 8;
pillar_height = 5;
num_hd_bins = 60;

% do not change, dependent on other variables/environmental geometry/default setting
wall_width = floor_width; pillar_width = floor_width/5;
viewbin_offset = 2;
num_place_bins = floor_width^2;
num_view_bins = viewbin_offset + 2*floor_width^2 + wall_height*wall_width*4 + 4*pillar_height*pillar_width*4;

% Code adapted from plotgridmap.m
exploded = false;  % exploded view, i.e. floor and ceiling are separated from walls/pillars

floor_x = repmat(0:floor_width, floor_width+1, 1);
floor_y = flipud(repmat([0:floor_width]', 1, floor_width+1));
floor_z = zeros(floor_width+1,floor_width+1);

ceiling_x = floor_x;
ceiling_y = floor_y;

walls_x = repmat([0.*ones(1,wall_width) 0:wall_width-1 wall_width.*ones(1,wall_width) wall_width:-1:0], wall_height+1, 1);
walls_y = repmat([0:wall_width-1 wall_width.*ones(1,wall_width) wall_width:-1:1 0.*ones(1,wall_width+1)], wall_height+1, 1);

P1_x = repmat([3*pillar_width.*ones(1,pillar_width) 3*pillar_width:4*pillar_width-1 4*pillar_width.*ones(1,pillar_width) 4*pillar_width:-1:3*pillar_width], pillar_height+1, 1);
P1_y = repmat([pillar_width:2*pillar_width-1 2*pillar_width.*ones(1,pillar_width) 2*pillar_width:-1:pillar_width+1 pillar_width.*ones(1,pillar_width+1)], pillar_height+1, 1);

P2_x = repmat([pillar_width.*ones(1,pillar_width) pillar_width:2*pillar_width-1 2*pillar_width.*ones(1,pillar_width) 2*pillar_width:-1:pillar_width], pillar_height+1, 1);
P2_y = P1_y;

P3_x = P1_x;
P3_y = repmat([3*pillar_width:4*pillar_width-1 4*pillar_width.*ones(1,pillar_width) 4*pillar_width:-1:3*pillar_width+1 3*pillar_width.*ones(1,pillar_width+1)], pillar_height+1, 1);

P4_x = P2_x;
P4_y = P3_y;

if exploded
    ceiling_z = 3*wall_height.*ones(floor_width+1,floor_width+1);
    walls_z = repmat([2*wall_height:-1:wall_height]', 1, wall_width*4 + 1);
    PX_z = repmat([wall_height+pillar_height:-1:wall_height]', 1, pillar_width*4 + 1);
else
    ceiling_z = wall_height.*ones(floor_width+1,floor_width+1);
    walls_z = repmat([wall_height:-1:0]', 1, wall_width*4 + 1);
    PX_z = repmat([pillar_height:-1:0]', 1, pillar_width*4 + 1);
end

floor_last_bin = floor_width^2 + viewbin_offset;
floor = flipud(reshape(viewbin_offset+1:floor_last_bin, floor_width, floor_width)');

% ceiling follows floor mapping, top down view
ceiling_last_bin = floor_last_bin + floor_width^2;
ceiling = flipud(reshape(floor_last_bin+1:ceiling_last_bin, floor_width, floor_width)');

% from top down, slit walls at bottom left corner, open outwards.
% start from row closest to ground, rightwards, then climb rows
walls_last_bin = ceiling_last_bin + 4*wall_width*wall_height;
walls = flipud(reshape(ceiling_last_bin+1:walls_last_bin, wall_width*4, wall_height)');

% BL - bottom left, and so on, from top view, same slicing as walls
% pillar width 8, height 5
P1_last_bin = walls_last_bin + 4*pillar_width*pillar_height;
P2_last_bin = P1_last_bin + 4*pillar_width*pillar_height;
P3_last_bin = P2_last_bin + 4*pillar_width*pillar_height;
P1_BR = flipud(reshape(walls_last_bin+1:P1_last_bin, pillar_width*4, pillar_height)');
P2_BL = flipud(reshape(P1_last_bin+1:P2_last_bin, pillar_width*4, pillar_height)');
P3_TR = flipud(reshape(P2_last_bin+1:P3_last_bin, pillar_width*4, pillar_height)');
P4_TL = flipud(reshape(P3_last_bin+1:num_view_bins, pillar_width*4, pillar_height)');

switch model
    case 'phv'
        params = cell2mat(params(:, 1)')';
        place_params = params(:, 1:num_place_bins);
        hd_params = params(:, num_place_bins+1:num_place_bins+num_hd_bins);
        view_params = params(:, num_place_bins+num_hd_bins+1:num_place_bins+num_hd_bins+num_view_bins);
        
        similarity_scores = similarity_scores{1};
        place_fit = similarity_scores(:,1);
        hd_fit = similarity_scores(:,2);
        view_fit = similarity_scores(:,3);
        
    case 'ph'
        params = cell2mat(params(:, 2)')';
        place_params = params(:, 1:num_place_bins);
        hd_params = params(:, num_place_bins+1:num_place_bins+num_hd_bins);
        
        similarity_scores = similarity_scores{2};
        place_fit = similarity_scores(:,1);
        hd_fit = similarity_scores(:,2);
        
    case 'pv'
        params = cell2mat(params(:, 3)')';
        place_params = params(:, 1:num_place_bins);
        view_params = params(:, num_place_bins+1:num_place_bins+num_view_bins);
        
        similarity_scores = similarity_scores{3};
        place_fit = similarity_scores(:,1);
        view_fit = similarity_scores(:,2);
        
    case 'hv'
        params = cell2mat(params(:, 4)')';
        hd_params = params(:, 1:num_hd_bins);
        view_params = params(:, num_hd_bins+1:num_hd_bins+num_view_bins);
        
        similarity_scores = similarity_scores{4};
        hd_fit = similarity_scores(:,1);
        view_fit = similarity_scores(:,2);
        
    case 'place'
        place_params = cell2mat(params(:, 5)')';
        place_fit = similarity_scores{5};
        
    case 'headdirection'
        hd_params = cell2mat(params(:, 6)')';
        hd_fit = similarity_scores{6};
        
    case 'spatialview'
        view_params = cell2mat(params(:, 7)')';
        view_fit = similarity_scores{7};
end

if strcmp(model, 'place') || strcmp(model, 'ph') || strcmp(model, 'pv') || strcmp(model, 'phv')
    fp = figure('Name','Place plot');
    axLims = zeros(num_folds, 2);
    
    for fc = 1:num_folds
        ratemap = exp(place_params(fc,:))/tbin_size;
        if filter_unoccupied
            ratemap(ratemap == exp(large_negative_number)/tbin_size) = NaN;
        end

        subplot(subplot_rows, subplot_cols, fc);
        surf(floor_x, floor_y, floor_z, flipud(reshape(ratemap(1:num_place_bins), floor_width, floor_width)'));
        alpha 1; shading flat;
        zlim([0,1]);
        view(-35,20);
        colormap jet;
        colorbar;
        
        rectangle('Position', [pillar_width, pillar_width, pillar_width, pillar_width], 'EdgeColor', 'k', 'LineWidth', 1);
        rectangle('Position', [pillar_width, 3*pillar_width, pillar_width, pillar_width], 'EdgeColor', 'k', 'LineWidth', 1);
        rectangle('Position', [3*pillar_width, pillar_width, pillar_width, pillar_width], 'EdgeColor', 'k', 'LineWidth', 1);
        rectangle('Position', [3*pillar_width, 3*pillar_width, pillar_width, pillar_width], 'EdgeColor', 'k', 'LineWidth', 1);
        
        text('String', ['ratemap fit: ' num2str(place_fit(fc))], 'Position', [-2, -2, -0.15], ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 11);
        axLims(fc, :) = [prctile(ratemap, pltRange(1)), prctile(ratemap, pltRange(2))];  % was previously caxis
    end
    
    if exist('axRange', 'var') && ~isempty(axRange)
        caxRange = axRange(1,:);
        axRange(1,:) = [];
    else
        caxRange = [0, max(axLims(:,2))];
    end
    for fc = 1:num_folds
        subplot(subplot_rows, subplot_cols, fc);
        caxis(caxRange);
    end
        annotation('textbox', [(subplot_cols-1)/subplot_cols, 0.25/subplot_rows, 0.2, 0.1], 'String', ['mean ratemap fit: ' num2str(nanmean(place_fit))], ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);
    if save
        saveas(fp, 'place_plot.fig');
    end
end

if strcmp(model, 'headdirection') || strcmp(model, 'ph') || strcmp(model, 'hv') || strcmp(model, 'phv')
    fh = figure('Name','Head direction plot');
    axLims = zeros(num_folds, 2);
    
    for fc = 1:num_folds
        ratemap = exp(hd_params(fc,:))/tbin_size;

        ax = subplot(subplot_rows, subplot_cols, fc);
        pax = polaraxes('Units', ax.Units, 'Position', ax.Position);
        polarplot(deg2rad((0:num_hd_bins)*360/num_hd_bins), [ratemap; ratemap(1)]);
        pax.ThetaZeroLocation = 'top';
        pax.ThetaDir = 'clockwise';
        set(ax, 'Visible', 'off');
        
        axLims(fc, :) = [prctile(ratemap, pltRange(1)), prctile(ratemap, pltRange(2))];  % was previously rlim
        text('String', ['ratemap fit: ' num2str(hd_fit(fc))], 'Position', [pi, 1.4*axLims(fc,2)], ...
            'HorizontalAlignment', 'center', 'FontSize', 11);
    end
    
    if exist('axRange', 'var') && ~isempty(axRange)
        caxRange = axRange(1,:);
        axRange(1,:) = [];
    else
        caxRange = [0, max(axLims(:,2))];
    end
    rlim(caxRange);
        annotation('textbox', [(subplot_cols-1)/subplot_cols, 0.25/subplot_rows, 0.2, 0.1], 'String', ['mean ratemap fit: ' num2str(nanmean(hd_fit))], ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);
    if save
        saveas(fh, 'hd_plot.fig');
    end
end

if strcmp(model, 'spatialview') || strcmp(model, 'pv') || strcmp(model, 'hv') || strcmp(model, 'phv')
    fv = figure('Name','View plot');
    axLims = zeros(num_folds, 2);
    
    for fc = 1:num_folds
        ratemap = exp(view_params(fc,:))/tbin_size;
        if filter_unoccupied
            ratemap(ratemap == exp(large_negative_number)/tbin_size) = NaN;
        end
        
        subplot(subplot_rows, subplot_cols, fc);
        
        % Plot floor
        surf(floor_x, floor_y, floor_z, flipud(reshape(ratemap(viewbin_offset+1:floor_last_bin), floor_width, floor_width)'));
        alpha 0.35; shading flat;
        hold on;

        % Plot ceiling and walls
        surf(ceiling_x, ceiling_y, ceiling_z, flipud(reshape(ratemap(floor_last_bin+1:ceiling_last_bin), floor_width, floor_width)'));
        alpha 0.35; shading flat;
        surf(walls_x, walls_y, walls_z, flipud(reshape(ratemap(ceiling_last_bin+1:walls_last_bin), wall_width*4, wall_height)'));      
        alpha 0.35; shading flat;

        % Plot pillars
        surf(P1_x, P1_y, PX_z, flipud(reshape(ratemap(walls_last_bin+1:P1_last_bin), pillar_width*4, pillar_height)'));
        alpha 0.35; shading flat;
        surf(P2_x, P2_y, PX_z, flipud(reshape(ratemap(P1_last_bin+1:P2_last_bin), pillar_width*4, pillar_height)'));
        alpha 0.35; shading flat;
        surf(P3_x, P3_y, PX_z, flipud(reshape(ratemap(P2_last_bin+1:P3_last_bin), pillar_width*4, pillar_height)'));
        alpha 0.35; shading flat;
        surf(P4_x, P4_y, PX_z, flipud(reshape(ratemap(P3_last_bin+1:num_view_bins), pillar_width*4, pillar_height)'));
        alpha 0.35; shading flat; 
        view(-35,20);
        colormap jet;
        colorbar;
        
        rectangle('Position', [pillar_width, pillar_width, pillar_width, pillar_width], 'EdgeColor', 'k', 'LineWidth', 1);
        rectangle('Position', [pillar_width, 3*pillar_width, pillar_width, pillar_width], 'EdgeColor', 'k', 'LineWidth', 1);
        rectangle('Position', [3*pillar_width, pillar_width, pillar_width, pillar_width], 'EdgeColor', 'k', 'LineWidth', 1);
        rectangle('Position', [3*pillar_width, 3*pillar_width, pillar_width, pillar_width], 'EdgeColor', 'k', 'LineWidth', 1);
        
        text('String', ['ratemap fit: ' num2str(view_fit(fc))], 'Position', [-2, -2, -1.5], ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 11);
        hold off;
        axLims(fc, :) = [prctile(ratemap, pltRange(1)), prctile(ratemap, pltRange(2))];  % was previously caxis
    end
    
    if exist('axRange', 'var') && ~isempty(axRange)
        caxRange = axRange(1,:);
        axRange(1,:) = [];
    else
        caxRange = [0, max(axLims(:,2))];
    end
    for fc = 1:num_folds
        subplot(subplot_rows, subplot_cols, fc);
        caxis(caxRange);
    end
    annotation('textbox', [(subplot_cols-1)/subplot_cols, 0.25/subplot_rows, 0.2, 0.1], 'String', ['mean ratemap fit: ' num2str(nanmean(view_fit))], ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);
    if save
        saveas(fv, 'view_plot.fig');
    end
end

end


function [rows, cols] = getSubplotGridSize(numSubplots)
    % Calculate the number of rows and columns for the subplot grid
    rows = ceil(sqrt(numSubplots));
    cols = ceil(numSubplots / rows);
end
