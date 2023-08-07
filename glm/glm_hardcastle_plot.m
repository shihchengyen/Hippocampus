function glm_hardcastle_plot(hc_results, model)
%	Plots fitted params onto each variable space,
%   analogous to the tuning curve of the cell.
%   Plots all folds in a single figure for each model variable.

%	PARAMETERS:
%	hc_results - struct output of glm_hardcastle
%	model - 'place' / 'headdirection' / 'view' / 
%           'ph' / 'pv' / 'hv' / 'phv'

% Code adapted from plotgridmap.m
params = hc_results.params_consol;
tbin_size = hc_results.tbin_size;
num_folds = size(params, 1);
[subplot_rows, subplot_cols] = getSubplotGridSize(num_folds);

floor_x = repmat(0:40, 41, 1);
floor_y = flipud(repmat([0:40]', 1, 41));
floor_z = zeros(41,41);

ceiling_x = floor_x;
ceiling_y = floor_y;
ceiling_z = 8.*ones(41,41);

walls_x = repmat([0.*ones(1,40) 0:39 40.*ones(1,40) 40:-1:0], 9, 1);
walls_y = repmat([0:39 40.*ones(1,40) 40:-1:1 0.*ones(1,41)], 9, 1);
walls_z = repmat([8:-1:0]', 1, 40*4 + 1);

P1_x = repmat([24.*ones(1,8) 24:31 32.*ones(1,8) 32:-1:24], 6, 1);
P1_y = repmat([8:15 16.*ones(1,8) 16:-1:9 8.*ones(1,9)], 6, 1);
PX_z = repmat([5:-1:0]', 1, 8*4 + 1);

P2_x = repmat([8.*ones(1,8) 8:15 16.*ones(1,8) 16:-1:8], 6, 1);
P2_y = P1_y;

P3_x = P1_x;
P3_y = repmat([24:31 32.*ones(1,8) 32:-1:25 24.*ones(1,9)], 6, 1);

P4_x = P2_x;
P4_y = P3_y;

floor = flipud(reshape(3:3+1600-1, 40, 40)');

% ceiling follows floor mapping, top down view
ceiling = flipud(reshape(1603:1603+1600-1, 40, 40)');

% from top down, slit walls at bottom left corner, open outwards.
% start from row closest to ground, rightwards, then climb rows
walls = flipud(reshape(3203:3203+1280-1, 40*4, 8)');

% BL - bottom left, and so on, from top view, same slicing as walls
% pillar width 8, height 5
P1_BR = flipud(reshape(4483:4483+160-1, 8*4, 5)');
P2_BL = flipud(reshape(4643:4643+160-1, 8*4, 5)');
P3_TR = flipud(reshape(4803:4803+160-1, 8*4, 5)');
P4_TL = flipud(reshape(4963:4963+160-1, 8*4, 5)');

switch model
    case 'phv'
        params = cell2mat(params(:, 1)')';
        place_params = params(:, 1:1600);
        hd_params = params(:, 1601:1600+60);
        view_params = params(:, 1601+60:1600+60+5122);
    case 'ph'
        params = cell2mat(params(:, 2)')';
        place_params = params(:, 1:1600);
        hd_params = params(:, 1601:1600+60);
    case 'pv'
        params = cell2mat(params(:, 3)')';
        place_params = params(:, 1:1600);
        view_params = params(:, 1601:1600+5122);
    case 'hv'
        params = cell2mat(params(:, 4)')';
        hd_params = params(:, 1:60);
        view_params = params(:, 61:60+5122);
    case 'place'
        place_params = cell2mat(params(:, 5)')';
    case 'headdirection'
        hd_params = cell2mat(params(:, 6)')';
    case 'view'
        view_params = cell2mat(params(:, 7)')';
end

if strcmp(model, 'place') || strcmp(model, 'ph') || strcmp(model, 'pv') || strcmp(model, 'phv')
    fp = figure('Name','Place plot');
    for fc = 1:num_folds
        ratemap = nan(1600,1);
        for k = 1:size(ratemap,1)
            ratemap(k) = exp(place_params(fc, k))/tbin_size;
        end
        
        subplot(subplot_rows, subplot_cols, fc);
        surf(floor_x, floor_y, floor_z, flipud(reshape(ratemap(1:1600), 40, 40)'));
        alpha 1; shading flat;
        view(-35,20);
        colormap jet;
        colorbar;
    end
    saveas(fp, 'place_plot.fig');
end

if strcmp(model, 'headdirection') || strcmp(model, 'ph') || strcmp(model, 'hv') || strcmp(model, 'phv')
    fh = figure('Name','Head direction plot');
    for fc = 1:num_folds
        ratemap = nan(60,1);
        for k = 1:size(ratemap,1)
            ratemap(k) = exp(hd_params(fc, k))/tbin_size;
        end

        ax = subplot(subplot_rows, subplot_cols, fc);
        pax = polaraxes('Units', ax.Units, 'Position', ax.Position);
        polarplot(deg2rad((0:60)*360/60), [ratemap; ratemap(1)]);
        pax.ThetaZeroLocation = 'top';
        pax.ThetaDir = 'clockwise';
        set(ax, 'Visible', 'off');
    end
    saveas(fh, 'hd_plot.fig');
end

if strcmp(model, 'view') || strcmp(model, 'pv') || strcmp(model, 'hv') || strcmp(model, 'phv')
    fv = figure('Name','View plot');
    for fc = 1:num_folds
        ratemap = nan(5122,1);
        for k = 1:size(ratemap,1)
            ratemap(k) = exp(view_params(fc, k))/tbin_size;
        end

        subplot(subplot_rows, subplot_cols, fc);
        
        % Plot floor
        surf(floor_x, floor_y, floor_z, flipud(reshape(ratemap(3:1600+3-1), 40, 40)'));
        alpha 0.35; shading flat;
        hold on;

        % Plot ceiling and walls
        surf(ceiling_x, ceiling_y, ceiling_z, flipud(reshape(ratemap(1603:1603+1600-1), 40, 40)'));
        alpha 0.35; shading flat;
        surf(walls_x, walls_y, walls_z, flipud(reshape(ratemap(3203:3203+1280-1), 40*4, 8)'));      
        alpha 0.35; shading flat;

        % Plot pillars
        surf(P1_x, P1_y, PX_z, flipud(reshape(ratemap(4483:4483+160-1), 8*4, 5)'));
        alpha 0.35; shading flat;
        surf(P2_x, P2_y, PX_z, flipud(reshape(ratemap(4643:4643+160-1), 8*4, 5)'));
        alpha 0.35; shading flat;
        surf(P3_x, P3_y, PX_z, flipud(reshape(ratemap(4803:4803+160-1), 8*4, 5)'));
        alpha 0.35; shading flat;
        surf(P4_x, P4_y, PX_z, flipud(reshape(ratemap(4963:4963+160-1), 8*4, 5)'));
        alpha 0.35; shading flat; 
        view(-35,20);
        colormap jet;
        colorbar;
    end
    saveas(fv, 'view_plot.fig');
end

end


function [rows, cols] = getSubplotGridSize(numSubplots)
    % Calculate the number of rows and columns for the subplot grid
    rows = ceil(sqrt(numSubplots));
    cols = ceil(numSubplots / rows);
end
