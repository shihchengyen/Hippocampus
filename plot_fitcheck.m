function plot_fitcheck()
% Creates an interactive GUI to browse through all cells with smoothed outlines
% Use arrow keys or buttons to navigate between cells

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
    fprintf('Got sample map from vpc in workspace\n');
elseif evalin('base', 'exist(''vmp'', ''var'')')
    vmp = evalin('base', 'vmp');
    sample_map = vmp.data.maps_adsm(1,:);
    fprintf('Got sample map from vmp in workspace\n');
else
    % Try loading from place_selective_cells.txt
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
                    fprintf('Got sample map from first cell in place_selective_cells.txt\n');
                elseif isfield(vpc_temp, 'vpc')
                    sample_map = vpc_temp.vpc.data.maps_adsm(1,:);
                    fprintf('Got sample map from first cell in place_selective_cells.txt\n');
                end
            end
        end
    end
    
    if isempty(sample_map)
        error('Cannot find sample map. Need vpc/vmp in workspace or place_selective_cells.txt');
    end
end

grid_size = 40;

% Extract cell data
cell_outlines = {};
for i = 1:size(data, 1)
    rc = data{i, 1};
    pi = data{i, 2};
    pv = data{i, 3};
    if ~isempty(rc) && iscell(rc) && ~isempty(pi) && ~isempty(pv)
        if size(data, 2) >= 5
            cell_outlines{end+1} = {rc, pi, pv, data{i,4}, data{i,5}};
        else
            cell_outlines{end+1} = {rc, pi, pv, NaN, false};
        end
    end
end

fprintf('Loaded %d cells for inspection\n', length(cell_outlines));

% Prepare maze background
map_2d = reshape(sample_map, grid_size, grid_size);
rgb_image = zeros(grid_size, grid_size, 3);
for r = 1:grid_size
    for c = 1:grid_size
        if isnan(map_2d(r, c))
            rgb_image(r, c, :) = [1 1 1];  % White
        else
            rgb_image(r, c, :) = [0.95 0.95 0.95];  % Light gray
        end
    end
end

% Load cell directories for labeling
cell_dirs = {};
if exist('place_selective_cells.txt', 'file')
    fid = fopen('place_selective_cells.txt', 'r');
    while ~feof(fid)
        line = fgetl(fid);
        if ischar(line) && ~isempty(line)
            cell_dirs{end+1} = line;
        end
    end
    fclose(fid);
end

% Create GUI
fig = figure('Position', [100 100 900 850], 'Name', 'Cell Inspector - Smoothed Outlines', ...
    'KeyPressFcn', @keyPress);

% Create axis for plot
ax = axes('Parent', fig, 'Position', [0.05 0.15 0.9 0.8]);

% Create navigation panel
uicontrol('Style', 'pushbutton', 'String', '<< Previous', ...
    'Units', 'normalized', 'Position', [0.1 0.05 0.15 0.05], ...
    'Callback', @previousCell, 'FontSize', 10, 'FontWeight', 'bold');

uicontrol('Style', 'pushbutton', 'String', 'Next >>', ...
    'Units', 'normalized', 'Position', [0.75 0.05 0.15 0.05], ...
    'Callback', @nextCell, 'FontSize', 10, 'FontWeight', 'bold');

% Cell counter text
cell_text = uicontrol('Style', 'text', 'String', '', ...
    'Units', 'normalized', 'Position', [0.35 0.05 0.3 0.05], ...
    'FontSize', 12, 'FontWeight', 'bold');

% Store data in figure
setappdata(fig, 'current_cell', 1);
setappdata(fig, 'total_cells', length(cell_outlines));
setappdata(fig, 'cell_outlines', cell_outlines);
setappdata(fig, 'rgb_image', rgb_image);
setappdata(fig, 'map_2d', map_2d);
setappdata(fig, 'grid_size', grid_size);
setappdata(fig, 'cell_dirs', cell_dirs);
setappdata(fig, 'ax', ax);
setappdata(fig, 'cell_text', cell_text);

% Plot first cell
plotCell(fig);

fprintf('\nGUI Controls:\n');
fprintf('  - Use << Previous / Next >> buttons to navigate\n');
fprintf('  - Use Left/Right arrow keys to navigate\n');
fprintf('  - Close window when done\n');

    function keyPress(src, event)
        if strcmp(event.Key, 'rightarrow')
            nextCell(src, event);
        elseif strcmp(event.Key, 'leftarrow')
            previousCell(src, event);
        end
    end

    function nextCell(~, ~)
        current = getappdata(fig, 'current_cell');
        total = getappdata(fig, 'total_cells');
        if current < total
            setappdata(fig, 'current_cell', current + 1);
            plotCell(fig);
        end
    end

    function previousCell(~, ~)
        current = getappdata(fig, 'current_cell');
        if current > 1
            setappdata(fig, 'current_cell', current - 1);
            plotCell(fig);
        end
    end

    function plotCell(fig)
        % Get data from figure
        current_cell = getappdata(fig, 'current_cell');
        total_cells = getappdata(fig, 'total_cells');
        cell_outlines = getappdata(fig, 'cell_outlines');
        rgb_image = getappdata(fig, 'rgb_image');
        map_2d = getappdata(fig, 'map_2d');
        grid_size = getappdata(fig, 'grid_size');
        cell_dirs = getappdata(fig, 'cell_dirs');
        ax = getappdata(fig, 'ax');
        cell_text = getappdata(fig, 'cell_text');
        
        % Clear axis
        cla(ax);
        
        % Plot background
        image(ax, rgb_image);
        set(ax, 'YDir', 'reverse');
        axis(ax, 'equal', 'tight');
        hold(ax, 'on');
        
        % Get current cell data
        region_cells = cell_outlines{current_cell}{1};
        
        cell_color = [0.2 0.4 0.8];  % Blue
        
        % Plot all regions for this cell
        for region_idx = 1:length(region_cells)
            current_region = region_cells{region_idx};
            [rows, cols] = ind2sub([grid_size, grid_size], current_region);
            
            if length(rows) < 3
                continue;
            end
            
            % Create mask for filled region
            mask = false(grid_size, grid_size);
            for i = 1:length(rows)
                mask(rows(i), cols(i)) = true;
            end
            
            % FILTER OUT pixels that are outside valid maze (in walls)
            valid_mask = false(grid_size, grid_size);
            for r = 1:grid_size
                for c = 1:grid_size
                    if mask(r, c) && ~isnan(map_2d(r, c))
                        valid_mask(r, c) = true;
                    end
                end
            end
            
            mask = valid_mask;
            
            % Fill the actual pixel region with transparency (0.3 alpha)
            [fill_rows, fill_cols] = find(mask);
            for i = 1:length(fill_rows)
                rectangle('Position', [fill_cols(i)-0.5, fill_rows(i)-0.5, 1, 1], ...
                    'FaceColor', [cell_color, 0.3], 'EdgeColor', 'none', 'Parent', ax);
            end
            
            % Create smoothed boundary outline
            boundaries = bwboundaries(mask, 'noholes');
            
            if isempty(boundaries)
                continue;
            end
            
            % Use the longest boundary (main outline)
            [~, max_idx] = max(cellfun(@(x) size(x,1), boundaries));
            boundary = boundaries{max_idx};
            
            boundary_y = boundary(:,1);  % rows
            boundary_x = boundary(:,2);  % cols
            
            % Ensure closure
            if boundary_x(end) ~= boundary_x(1) || boundary_y(end) ~= boundary_y(1)
                boundary_x = [boundary_x; boundary_x(1)];
                boundary_y = [boundary_y; boundary_y(1)];
            end
            
            % Apply smoothing
            window_size = 5;  % Adjust for different smoothness
            
            if length(boundary_x) > window_size
                % Pad the arrays for wraparound smoothing
                pad_size = ceil(window_size/2);
                padded_x = [boundary_x(end-pad_size:end-1); boundary_x; boundary_x(2:pad_size+1)];
                padded_y = [boundary_y(end-pad_size:end-1); boundary_y; boundary_y(2:pad_size+1)];
                
                % Smooth with padding
                smooth_x_full = movmean(padded_x, window_size);
                smooth_y_full = movmean(padded_y, window_size);
                
                % Extract the middle portion
                smooth_x = smooth_x_full(pad_size+1:end-pad_size);
                smooth_y = smooth_y_full(pad_size+1:end-pad_size);
                
                % Ensure loop closure
                smooth_x(end) = smooth_x(1);
                smooth_y(end) = smooth_y(1);
                
                % Plot the smoothed boundary
                plot(ax, smooth_x, smooth_y, '-', 'Color', cell_color, 'LineWidth', 2.5);
            else
                % If region too small, just plot the boundary
                plot(ax, boundary_x, boundary_y, '-', 'Color', cell_color, 'LineWidth', 2.5);
            end
        end
        
        hold(ax, 'off');
        
        % Create title with session info
        title_text = sprintf('Cell %d of %d - %d Place Fields', ...
            current_cell, total_cells, length(region_cells));
        
        % Add session directory info if available
        if current_cell <= length(cell_dirs)
            cell_path = cell_dirs{current_cell};
            
            % Extract date
            date_match = regexp(cell_path, '\d{8}', 'match');
            session_date = '';
            if ~isempty(date_match)
                session_date = date_match{1};
            end
            
            % Extract session, array, channel, cell numbers
            session_match = regexp(cell_path, 'session(\d+)', 'tokens');
            array_match = regexp(cell_path, 'array(\d+)', 'tokens');
            channel_match = regexp(cell_path, 'channel(\d+)', 'tokens');
            cell_match = regexp(cell_path, 'cell(\d+)', 'tokens');
            
            title_parts = {};
            
            if ~isempty(session_date)
                title_parts{end+1} = session_date;
            end
            
            if ~isempty(session_match)
                title_parts{end+1} = sprintf('session %s', session_match{1}{1});
            end
            
            if ~isempty(array_match)
                title_parts{end+1} = sprintf('array %s', array_match{1}{1});
            end
            
            if ~isempty(channel_match)
                title_parts{end+1} = sprintf('channel %s', channel_match{1}{1});
            end
            
            if ~isempty(cell_match)
                title_parts{end+1} = sprintf('cell %s', cell_match{1}{1});
            end
            
            if ~isempty(title_parts)
                session_info = strjoin(title_parts, ' | ');
                title_text = sprintf('Cell %d of %d: %s\n%d Place Fields', ...
                    current_cell, total_cells, session_info, length(region_cells));
            end
        end
        
        title(ax, title_text, 'FontSize', 10, 'FontWeight', 'bold', 'Interpreter', 'none');
        
        % Update counter text
        set(cell_text, 'String', sprintf('Cell %d / %d', current_cell, total_cells));
    end

end