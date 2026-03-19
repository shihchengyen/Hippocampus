function plot_3d_overlay(varargin)
% PLOT_3D_OVERLAY Creates 3D visualization of place fields stacked vertically
% Each cell gets its own z-level, creating a 3D overlay effect
%
% USAGE:
%   plot_3d_overlay()                    - 3D stack of ALL cells (default)
%   plot_3d_overlay('pixel_perfect')     - Use pixel-perfect edges

% Parse input arguments
use_smoothed_outline = true;

if nargin > 0
    for i = 1:nargin
        if strcmpi(varargin{i}, 'pixel_perfect')
            use_smoothed_outline = false;
        end
    end
end

% Check if data variable exists
if ~evalin('base', 'exist(''data'', ''var'')')
    error('No ''data'' variable found in workspace. Run ProcessDirs first.');
end

% Get data from base workspace
data = evalin('base', 'data');

% Get sample map for maze structure
sample_map = [];
if evalin('base', 'exist(''vpc'', ''var'')')
    vpc = evalin('base', 'vpc');
    sample_map = vpc.data.maps_adsm(1,:);
elseif evalin('base', 'exist(''vmp'', ''var'')')
    vmp = evalin('base', 'vmp');
    sample_map = vmp.data.maps_adsm(1,:);
else
    if exist('place_selective_cells.txt', 'file')
        fid = fopen('place_selective_cells.txt', 'r');
        cell_dirs = {};
        while ~feof(fid)
            line = fgetl(fid);
            if ischar(line) && ~isempty(line)
                cell_dirs{end+1} = line;
            end
        end
        fclose(fid);
        
        if ~isempty(cell_dirs)
            first_cell_dir = cell_dirs{1};
            vmpc_file = fullfile(first_cell_dir,'vmpc.mat');
            if exist(vmpc_file, 'file')
                vpc_temp = load(vmpc_file);
                if isfield(vpc_temp, 'vmp')
                    sample_map = vpc_temp.vmp.data.maps_adsm(1,:);
                elseif isfield(vpc_temp, 'vpc')
                    sample_map = vpc_temp.vpc.data.maps_adsm(1,:);
                end
            end
        end
    end
end

if isempty(sample_map)
    error('Cannot find sample map.');
end

grid_size = 40;

% Extract cell data - just get all cells
cell_outlines = {};
for i = 1:size(data, 1)
    rc = data{i, 1};
    pi = data{i, 2};
    pv = data{i, 3};
    if ~isempty(rc) && iscell(rc) && ~isempty(pi) && ~isempty(pv)
        cell_outlines{end+1} = {rc, pi, pv};
    end
end

fprintf('Creating 3D visualization with ALL %d cells\n', length(cell_outlines));

if isempty(cell_outlines)
    warning('No cells to plot.');
    return;
end

% Use all cells
cells_to_plot = cell_outlines;

% Generate colors
colors = hsv(length(cells_to_plot));

% Prepare maze background
map_2d = reshape(sample_map, grid_size, grid_size);

%% CREATE 3D VISUALIZATION
figure('Position', [100 100 1000 900], 'Name', '3D Place Field Overlay');

% Set up 3D axes
ax = axes('Position', [0.05 0.05 0.9 0.85]);
hold(ax, 'on');
grid(ax, 'on');
view(ax, 3);

% Set z-spacing between layers
z_spacing = 1.5;

% Draw maze base at z=0
for r = 1:grid_size
    for c = 1:grid_size
        if ~isnan(map_2d(r, c))
            % Draw maze floor as light gray patches
            patch(ax, [c-0.5 c+0.5 c+0.5 c-0.5], [r-0.5 r-0.5 r+0.5 r+0.5], ...
                [0 0 0 0], [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        end
    end
end

% Plot each cell at different z-level
for cell_idx = 1:length(cells_to_plot)
    region_cells = cells_to_plot{cell_idx}{1};
    
    if ~iscell(region_cells) || isempty(region_cells)
        continue;
    end
    
    cell_color = colors(cell_idx, :);
    z_level = cell_idx * z_spacing;
    
    % Plot all regions for this cell at the same z-level
    for region_idx = 1:length(region_cells)
        current_region = region_cells{region_idx};
        [rows, cols] = ind2sub([grid_size, grid_size], current_region);
        
        if length(rows) < 3
            continue;
        end
        
        % Create mask
        mask = false(grid_size, grid_size);
        for i = 1:length(rows)
            mask(rows(i), cols(i)) = true;
        end
        
        % Filter valid maze pixels
        valid_mask = false(grid_size, grid_size);
        for r = 1:grid_size
            for c = 1:grid_size
                if mask(r, c) && ~isnan(map_2d(r, c))
                    valid_mask(r, c) = true;
                end
            end
        end
        mask = valid_mask;
        
        % Draw filled region at this z-level
        [fill_rows, fill_cols] = find(mask);
        for i = 1:length(fill_rows)
            r = fill_rows(i);
            c = fill_cols(i);
            patch(ax, [c-0.5 c+0.5 c+0.5 c-0.5], [r-0.5 r-0.5 r+0.5 r+0.5], ...
                [z_level z_level z_level z_level], cell_color, ...
                'EdgeColor', 'none', 'FaceAlpha', 0.6);
        end
        
        if use_smoothed_outline
            % SMOOTHED OUTLINE in 3D
            boundaries = bwboundaries(mask, 'noholes');
            
            if ~isempty(boundaries)
                [~, max_idx] = max(cellfun(@(x) size(x,1), boundaries));
                boundary = boundaries{max_idx};
                
                boundary_y = boundary(:,1);
                boundary_x = boundary(:,2);
                
                if boundary_x(end) ~= boundary_x(1) || boundary_y(end) ~= boundary_y(1)
                    boundary_x = [boundary_x; boundary_x(1)];
                    boundary_y = [boundary_y; boundary_y(1)];
                end
                
                window_size = 5;
                
                if length(boundary_x) > window_size
                    pad_size = ceil(window_size/2);
                    padded_x = [boundary_x(end-pad_size:end-1); boundary_x; boundary_x(2:pad_size+1)];
                    padded_y = [boundary_y(end-pad_size:end-1); boundary_y; boundary_y(2:pad_size+1)];
                    
                    smooth_x_full = movmean(padded_x, window_size);
                    smooth_y_full = movmean(padded_y, window_size);
                    
                    smooth_x = smooth_x_full(pad_size+1:end-pad_size);
                    smooth_y = smooth_y_full(pad_size+1:end-pad_size);
                    
                    smooth_x(end) = smooth_x(1);
                    smooth_y(end) = smooth_y(1);
                    
                    z_coords = ones(size(smooth_x)) * z_level;
                    plot3(ax, smooth_x, smooth_y, z_coords, '-', 'Color', cell_color, 'LineWidth', 2);
                else
                    z_coords = ones(size(boundary_x)) * z_level;
                    plot3(ax, boundary_x, boundary_y, z_coords, '-', 'Color', cell_color, 'LineWidth', 2);
                end
            end
        else
            % PIXEL-PERFECT OUTLINE in 3D
            [fill_rows, fill_cols] = find(mask);
            
            for i = 1:length(fill_rows)
                r = fill_rows(i);
                c = fill_cols(i);
                
                if r == 1 || ~mask(r-1, c)
                    plot3(ax, [c-0.5, c+0.5], [r-0.5, r-0.5], [z_level, z_level], ...
                        '-', 'Color', cell_color, 'LineWidth', 2);
                end
                
                if r == grid_size || ~mask(r+1, c)
                    plot3(ax, [c-0.5, c+0.5], [r+0.5, r+0.5], [z_level, z_level], ...
                        '-', 'Color', cell_color, 'LineWidth', 2);
                end
                
                if c == 1 || ~mask(r, c-1)
                    plot3(ax, [c-0.5, c-0.5], [r-0.5, r+0.5], [z_level, z_level], ...
                        '-', 'Color', cell_color, 'LineWidth', 2);
                end
                
                if c == grid_size || ~mask(r, c+1)
                    plot3(ax, [c+0.5, c+0.5], [r-0.5, r+0.5], [z_level, z_level], ...
                        '-', 'Color', cell_color, 'LineWidth', 2);
                end
            end
        end
    end
end

hold(ax, 'off');

% Set axes properties
xlabel(ax, 'X Position', 'FontSize', 12, 'FontWeight', 'bold');
ylabel(ax, 'Y Position', 'FontSize', 12, 'FontWeight', 'bold');
zlabel(ax, 'Cell Layer', 'FontSize', 12, 'FontWeight', 'bold');

xlim(ax, [0 grid_size+1]);
ylim(ax, [0 grid_size+1]);
zlim(ax, [0 (length(cells_to_plot)+1)*z_spacing]);

title(ax, sprintf('3D Place Field Overlay: %d Cells', length(cells_to_plot)), ...
    'FontSize', 14, 'FontWeight', 'bold');
filename = 'cells_3d_overlay_all.png';

% Add lighting for better 3D effect
camlight('headlight');
lighting gouraud;

% Set nice viewing angle
view(ax, 45, 30);

saveas(gcf, filename);
fprintf('Saved to %s\n', filename);
fprintf('\nTip: Rotate the view using the mouse to see different angles!\n');

end