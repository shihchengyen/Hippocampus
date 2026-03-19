function plot_outline_overlay(varargin)
% PLOT_OUTLINE_OVERLAY Creates plots from data variable in workspace
% SEQUENCE:
% 1. Identify significant place cells using individual SIC thresholds
% 2. Generate SIC histogram with thresholds (OPTIONAL - use 'plot_sic' flag)
% 3. Plot outlines for cells (significant by default, or all with 'plot_all' flag)
% 4. Generate histogram of place field counts per cell
%
% USAGE:
%   plot_outline_overlay()                    - Smoothed outlines, significant cells only
%   plot_outline_overlay('plot_sic')          - Also generate SIC histogram
%   plot_outline_overlay('pixel_perfect')     - Use pixel-perfect edges instead of smoothed
%   plot_outline_overlay('plot_all')          - Plot ALL cells instead of just significant
%   plot_outline_overlay('plot_sic', 'plot_all', 'pixel_perfect') - Combine options

% Parse input arguments
plot_sic_histogram = false;
use_smoothed_outline = true;  % Default to smoothed
plot_all_cells = false;        % Default to significant cells only

if nargin > 0
    for i = 1:nargin
        if strcmpi(varargin{i}, 'plot_sic')
            plot_sic_histogram = true;
        elseif strcmpi(varargin{i}, 'pixel_perfect')
            use_smoothed_outline = false;
        elseif strcmpi(varargin{i}, 'plot_all')
            plot_all_cells = true;
        end
    end
end

% Check if data variable exists
if ~evalin('base', 'exist(''data'', ''var'')')
    error('No ''data'' variable found in workspace. Run ProcessDirs first.');
end

% Get data from base workspace
data = evalin('base', 'data');

% Validate data structure
if isempty(data) || size(data, 2) < 3
    error('Data variable has incorrect structure. Expected Nx3, Nx4, or Nx5 cell array.');
end

% Check if SIC values and significance are included
has_sic = (size(data, 2) >= 4);
has_significance = (size(data, 2) >= 5);

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
    
    % If still no sample map, try sessions
    if isempty(sample_map) && evalin('base', 'exist(''sessions'', ''var'')')
        sessions = evalin('base', 'sessions');
        if isfield(sessions, 'SessionDirs') && ~isempty(sessions.SessionDirs)
            if iscell(sessions.SessionDirs)
                first_cell_dir = sessions.SessionDirs{1};
            else
                first_cell_dir = sessions.SessionDirs;
            end
            vmpc_file = fullfile(first_cell_dir,'vmpc.mat');
            if exist(vmpc_file, 'file')
                vpc_temp = load(vmpc_file);
                if isfield(vpc_temp, 'vmp')
                    sample_map = vpc_temp.vmp.data.maps_adsm(1,:);
                    fprintf('Got sample map from sessions.SessionDirs\n');
                elseif isfield(vpc_temp, 'vpc')
                    sample_map = vpc_temp.vpc.data.maps_adsm(1,:);
                    fprintf('Got sample map from sessions.SessionDirs\n');
                end
            end
        end
    end
    
    if isempty(sample_map)
        error('Cannot find sample map. Need vpc/vmp in workspace, place_selective_cells.txt, or sessions.SessionDirs');
    end
end

grid_size = 40;

% Extract outlines from data structure
cell_outlines = {};
valid_cells = 0;

for i = 1:size(data, 1)
    rc = data{i, 1};  % region_cells
    pi = data{i, 2};  % peak_indices
    pv = data{i, 3};  % peak_values
    
    if ~isempty(rc) && iscell(rc) && ~isempty(pi) && ~isempty(pv)
        valid_cells = valid_cells + 1;
        if has_significance
            sic = data{i, 4};  % SIC value
            sig = data{i, 5};  % is_significant
            cell_outlines{valid_cells} = {rc, pi, pv, sic, sig};
        elseif has_sic
            sic = data{i, 4};  % SIC value
            cell_outlines{valid_cells} = {rc, pi, pv, sic, false};
        else
            cell_outlines{valid_cells} = {rc, pi, pv, NaN, false};
        end
    else
        fprintf('Skipping cell %d: empty or invalid data\n', i);
    end
end

fprintf('Found %d cells with valid outline data\n', valid_cells);

if valid_cells == 0
    warning('No cells with valid outline data found.');
    return;
end

%% STEP 1: DETERMINE SIGNIFICANCE USING INDIVIDUAL THRESHOLDS
fprintf('\n=== STEP 1: Identifying Significant Place Cells ===\n');

sic_values = zeros(length(cell_outlines), 1);
individual_thresholds = zeros(length(cell_outlines), 1);
is_significant = false(length(cell_outlines), 1);

% Extract SIC values
for i = 1:length(cell_outlines)
    sic_values(i) = cell_outlines{i}{4};
end

% Load individual thresholds
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
    fprintf('Loaded %d cell directories from place_selective_cells.txt\n', length(cell_dirs));
    
    % Load thresholds for each cell
    for i = 1:length(cell_outlines)
        if i > length(cell_dirs)
            fprintf('  Cell %d: No corresponding directory\n', i);
            continue;
        end
        
        cell_dir = cell_dirs{i};
        
        % Try multiple possible locations for vmpc.mat
        possible_paths = {
            fullfile(cell_dir, 'FiltAll', '1px', 'vmpc.mat'),
            fullfile(cell_dir, 'FiltAll', 'vmpc.mat'),
            fullfile(cell_dir, 'vmpc.mat')
        };
        
        vpc_cell = [];
        for p = 1:length(possible_paths)
            if exist(possible_paths{p}, 'file')
                try
                    vpc_temp = load(possible_paths{p});
                    if isfield(vpc_temp, 'vmp')
                        vpc_cell = vpc_temp.vmp;
                        break;
                    elseif isfield(vpc_temp, 'vpc')
                        vpc_cell = vpc_temp.vpc;
                        break;
                    end
                catch
                    continue;
                end
            end
        end
        
        if isempty(vpc_cell)
            fprintf('  Cell %d: Cannot find valid vmpc.mat\n', i);
            continue;
        end
        
        % Get individual threshold from shuffle data
        threshold_found = false;
        
        if isfield(vpc_cell.data, 'critsh_sm') && ~isempty(vpc_cell.data.critsh_sm)
            individual_thresholds(i) = prctile(vpc_cell.data.critsh_sm, 95);
            threshold_found = true;
        elseif isfield(vpc_cell.data, 'SICsh') && ~isempty(vpc_cell.data.SICsh)
            individual_thresholds(i) = prctile(vpc_cell.data.SICsh, 95);
            threshold_found = true;
        end
        
        % Determine significance
        if threshold_found && ~isnan(sic_values(i)) && individual_thresholds(i) > 0
            is_significant(i) = (sic_values(i) > individual_thresholds(i));
            if is_significant(i)
                fprintf('  Cell %d: ? SIGNIFICANT (SIC %.4f > threshold %.4f)\n', ...
                    i, sic_values(i), individual_thresholds(i));
            end
        end
    end
end

num_significant = sum(is_significant);
fprintf('\nIdentified %d significant place cells (%.1f%%)\n', ...
    num_significant, 100*num_significant/length(cell_outlines));

%% STEP 2: GENERATE SIC HISTOGRAM (OPTIONAL)
if plot_sic_histogram
    fprintf('\n=== STEP 2: Creating SIC Histogram ===\n');
    figure('Position', [50 100 1000 500], 'Name', 'SIC Histogram - All Cells');
    
    num_valid = sum(~isnan(sic_values));
    
    if num_valid > 0
        x_pos = 1:length(sic_values);
        
        % Plot SIC bars
        bar_handle = bar(x_pos, sic_values, 'FaceColor', 'flat');
        hold on;
        
        % Plot individual thresholds as red line with circles
        valid_idx = individual_thresholds > 0;
        if any(valid_idx)
            plot(x_pos(valid_idx), individual_thresholds(valid_idx), 'r-', 'LineWidth', 2.5);
            hScatter = scatter(x_pos(valid_idx), individual_thresholds(valid_idx), 80, 'r', 'filled', ...
                    'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        end
        
        % Add asterisks above significant cells
        sig_idx = is_significant & ~isnan(sic_values);
        if any(sig_idx)
            y_asterisk = zeros(sum(sig_idx), 1);
            sig_positions = x_pos(sig_idx);
            for i = 1:sum(sig_idx)
                idx = find(sig_idx, i);
                idx = idx(end);
                y_asterisk(i) = max(sic_values(idx), individual_thresholds(idx)) + 0.01;
            end
            hSig = plot(sig_positions, y_asterisk, 'g*', 'MarkerSize', 10, 'LineWidth', 2);
        end
        
        hold off;
        
        xlabel('Cell Number', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('SIC (bits/spike)', 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('SIC Values vs Individual 95%% Thresholds (ALL %d Cells)', length(sic_values)), ...
            'FontSize', 14, 'FontWeight', 'bold');
        grid on;
        xlim([0 length(sic_values)+1]);
        
        % Auto-scale y-axis
        valid_sic = sic_values(~isnan(sic_values));
        valid_thresh = individual_thresholds(individual_thresholds > 0);
        if ~isempty(valid_sic) || ~isempty(valid_thresh)
            all_values = [valid_sic; valid_thresh];
            y_min = min(all_values);
            y_max = max(all_values);
            y_range = y_max - y_min;
            ylim([max(0, y_min - 0.1*y_range), y_max + 0.15*y_range]);
        end
        
        % Add legend
        if exist('hScatter', 'var') && exist('hSig', 'var')
            legend([hScatter, hSig], {'Threshold', 'Significant'}, 'Location', 'best', 'FontSize', 11);
        end
        
        % Statistics text
        stats_text = sprintf(['Total: %d cells\nSignificant: %d (%.1f%%)\n' ...
            'Mean SIC: %.4f\nMedian SIC: %.4f'], ...
            length(sic_values), num_significant, 100*num_significant/length(sic_values), ...
            mean(valid_sic), median(valid_sic));
        
        annotation('textbox', [0.15, 0.75, 0.25, 0.15], ...
            'String', stats_text, ...
            'FitBoxToText', 'on', ...
            'BackgroundColor', 'white', ...
            'EdgeColor', 'black', ...
            'LineWidth', 1.5, ...
            'FontSize', 10, ...
            'FontWeight', 'bold');
        
        saveas(gcf, 'cells_sic_histogram_all_cells.png');
        fprintf('Saved to cells_sic_histogram_all_cells.png\n');
    end
else
    fprintf('\n=== STEP 2: SIC Histogram skipped (use plot_outline_overlay(''plot_sic'') to generate) ===\n');
end

%% STEP 3: FILTER TO SIGNIFICANT CELLS ONLY (OR USE ALL)
if plot_all_cells
    fprintf('\n=== STEP 3: Using ALL Cells for Plotting ===\n');
    cells_to_plot = cell_outlines;
    cell_indices_to_plot = 1:length(cell_outlines);
    fprintf('Processing all %d cells for plotting\n', length(cells_to_plot));
else
    fprintf('\n=== STEP 3: Filtering to Significant Cells Only ===\n');
    cells_to_plot = {};
    cell_indices_to_plot = [];
    for i = 1:length(cell_outlines)
        if is_significant(i)
            cells_to_plot{end+1} = cell_outlines{i};
            cell_indices_to_plot(end+1) = i;
        end
    end
    fprintf('Processing %d significant cells for plotting\n', length(cells_to_plot));
end

if isempty(cells_to_plot)
    warning('No cells found for plotting.');
    return;
end

% Generate colors for cells
colors = hsv(length(cells_to_plot));

% Prepare maze structure image
map_2d = reshape(sample_map, grid_size, grid_size);
rgb_image = zeros(grid_size, grid_size, 3);
for r = 1:grid_size
    for c = 1:grid_size
        if isnan(map_2d(r, c))
            rgb_image(r, c, :) = [1 1 1];  % White background
        else
            rgb_image(r, c, :) = [0.95 0.95 0.95];  % Light gray maze
        end
    end
end

%% STEP 4: PLOT ALL REGIONS
fprintf('\n=== STEP 4: Creating All Regions Plot ===\n');
if plot_all_cells
    fig_title = sprintf('All Regions: %d Cells (All)', length(cells_to_plot));
    filename = 'cells_outline_all_regions_all.png';
else
    fig_title = sprintf('All Regions: %d Cells (Significant)', length(cells_to_plot));
    filename = 'cells_outline_all_regions_significant.png';
end

figure('Position', [100 150 800 800], 'Name', 'All Regions');
image(rgb_image);
set(gca, 'YDir', 'reverse');
axis equal tight;
hold on;

for cell_idx = 1:length(cells_to_plot)
    region_cells = cells_to_plot{cell_idx}{1};
    
    if ~iscell(region_cells) || isempty(region_cells)
        continue;
    end
    
    cell_color = colors(cell_idx, :);
    
    for region_idx = 1:length(region_cells)
        current_region = region_cells{region_idx};
        [rows, cols] = ind2sub([grid_size, grid_size], current_region);
        
        if length(rows) < 3
            continue;
        end
        
        % Create mask for this region
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
        
        if use_smoothed_outline
            % SMOOTHED OUTLINE
            boundaries = bwboundaries(mask, 'noholes');
            
            if isempty(boundaries)
                continue;
            end
            
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
                
                plot(smooth_x, smooth_y, '-', 'Color', cell_color, 'LineWidth', 2.5);
            else
                plot(boundary_x, boundary_y, '-', 'Color', cell_color, 'LineWidth', 2.5);
            end
        else
            % PIXEL-PERFECT OUTLINE
            [fill_rows, fill_cols] = find(mask);
            
            for i = 1:length(fill_rows)
                r = fill_rows(i);
                c = fill_cols(i);
                
                if r == 1 || ~mask(r-1, c)
                    plot([c-0.5, c+0.5], [r-0.5, r-0.5], '-', 'Color', cell_color, 'LineWidth', 2.5);
                end
                
                if r == grid_size || ~mask(r+1, c)
                    plot([c-0.5, c+0.5], [r+0.5, r+0.5], '-', 'Color', cell_color, 'LineWidth', 2.5);
                end
                
                if c == 1 || ~mask(r, c-1)
                    plot([c-0.5, c-0.5], [r-0.5, r+0.5], '-', 'Color', cell_color, 'LineWidth', 2.5);
                end
                
                if c == grid_size || ~mask(r, c+1)
                    plot([c+0.5, c+0.5], [r-0.5, r+0.5], '-', 'Color', cell_color, 'LineWidth', 2.5);
                end
            end
        end
    end
end

hold off;
title(fig_title, 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, filename);
fprintf('Saved to %s\n', filename);

%% STEP 5: PLOT STRONGEST REGION ONLY
fprintf('\n=== STEP 5: Creating Strongest Region Plot ===\n');
if plot_all_cells
    fig_title = sprintf('Strongest Region Only: %d Cells (All)', length(cells_to_plot));
    filename = 'cells_outline_strongest_only_all.png';
else
    fig_title = sprintf('Strongest Region Only: %d Cells (Significant)', length(cells_to_plot));
    filename = 'cells_outline_strongest_only_significant.png';
end

figure('Position', [150 200 800 800], 'Name', 'Strongest Region');
image(rgb_image);
set(gca, 'YDir', 'reverse');
axis equal tight;
hold on;

for cell_idx = 1:length(cells_to_plot)
    region_cells = cells_to_plot{cell_idx}{1};
    peak_values = cells_to_plot{cell_idx}{3};
    
    if ~iscell(region_cells) || isempty(region_cells)
        continue;
    end
    
    % Find strongest region
    [~, strongest_idx] = max(peak_values);
    strongest_region = region_cells{strongest_idx};
    
    cell_color = colors(cell_idx, :);
    [rows, cols] = ind2sub([grid_size, grid_size], strongest_region);
    
    if length(rows) < 3
        continue;
    end
    
    % Create mask for this region
    mask = false(grid_size, grid_size);
    mask(sub2ind([grid_size, grid_size], rows, cols)) = true;
    
    % FILTER OUT pixels that are outside valid maze
    valid_mask = false(grid_size, grid_size);
    for r = 1:grid_size
        for c = 1:grid_size
            if mask(r, c) && ~isnan(map_2d(r, c))
                valid_mask(r, c) = true;
            end
        end
    end
    
    mask = valid_mask;
    
    if use_smoothed_outline
        % SMOOTHED OUTLINE
        boundaries = bwboundaries(mask, 'noholes');
        
        if isempty(boundaries)
            continue;
        end
        
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
            
            plot(smooth_x, smooth_y, '-', 'Color', cell_color, 'LineWidth', 2.5);
        else
            plot(boundary_x, boundary_y, '-', 'Color', cell_color, 'LineWidth', 2.5);
        end
    else
        % PIXEL-PERFECT OUTLINE
        [fill_rows, fill_cols] = find(mask);
        
        for i = 1:length(fill_rows)
            r = fill_rows(i);
            c = fill_cols(i);
            
            if r == 1 || ~mask(r-1, c)
                plot([c-0.5, c+0.5], [r-0.5, r-0.5], '-', 'Color', cell_color, 'LineWidth', 2.5);
            end
            
            if r == grid_size || ~mask(r+1, c)
                plot([c-0.5, c+0.5], [r+0.5, r+0.5], '-', 'Color', cell_color, 'LineWidth', 2.5);
            end
            
            if c == 1 || ~mask(r, c-1)
                plot([c-0.5, c-0.5], [r-0.5, r+0.5], '-', 'Color', cell_color, 'LineWidth', 2.5);
            end
            
            if c == grid_size || ~mask(r, c+1)
                plot([c+0.5, c+0.5], [r-0.5, r+0.5], '-', 'Color', cell_color, 'LineWidth', 2.5);
            end
        end
    end
end

hold off;
title(fig_title, 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, filename);
fprintf('Saved to %s\n', filename);

%% STEP 6: GENERATE PLACE FIELD COUNT HISTOGRAM
fprintf('\n=== STEP 6: Generating Place Field Count Histogram ===\n');

% Count place fields for each cell being plotted
num_place_fields = zeros(length(cells_to_plot), 1);
for i = 1:length(cells_to_plot)
    region_cells = cells_to_plot{i}{1};
    num_place_fields(i) = length(region_cells);
end

% Create bar plot (one bar per cell)
figure('Position', [200 250 1000 600], 'Name', 'Place Field Count per Cell');
x_pos = 1:length(num_place_fields);
bar(x_pos, num_place_fields, 'FaceColor', [0.3 0.6 0.9], 'EdgeColor', 'black', 'LineWidth', 1);

xlabel('Cell Number', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Number of Place Fields', 'FontSize', 12, 'FontWeight', 'bold');
if plot_all_cells
    title(sprintf('Place Field Count per Cell (n=%d cells, ALL)', length(cells_to_plot)), ...
        'FontSize', 14, 'FontWeight', 'bold');
else
    title(sprintf('Place Field Count per Cell (n=%d cells, Significant)', length(cells_to_plot)), ...
        'FontSize', 14, 'FontWeight', 'bold');
end
grid on;

% Set axis ranges with integer ticks only
xlim([0 length(num_place_fields)+1]);
y_max = max(num_place_fields);
ylim([0 y_max + 1]);
set(gca, 'YTick', 0:1:y_max+1);
set(gca, 'XTick', 1:1:length(num_place_fields));

% Add statistics box
mean_fields = mean(num_place_fields);
median_fields = median(num_place_fields);
stats_text = sprintf(['Total cells: %d\nTotal fields: %d\n' ...
    'Mean: %.2f\nMedian: %.0f\nMode: %d\nRange: %d-%d'], ...
    length(num_place_fields), sum(num_place_fields), ...
    mean_fields, median_fields, mode(num_place_fields), ...
    min(num_place_fields), max(num_place_fields));

annotation('textbox', [0.70, 0.65, 0.25, 0.25], ...
    'String', stats_text, ...
    'FitBoxToText', 'on', ...
    'BackgroundColor', 'white', ...
    'EdgeColor', 'black', ...
    'LineWidth', 1.5, ...
    'FontSize', 10, ...
    'FontWeight', 'bold');

if plot_all_cells
    saveas(gcf, 'place_field_count_histogram_all.png');
    fprintf('Saved to place_field_count_histogram_all.png\n');
else
    saveas(gcf, 'place_field_count_histogram.png');
    fprintf('Saved to place_field_count_histogram.png\n');
end

%% FINAL SUMMARY
fprintf('\n=== OVERALL SUMMARY ===\n');
fprintf('Total cells analyzed: %d\n', length(cell_outlines));
fprintf('Significant place cells: %d (%.1f%%)\n', num_significant, 100*num_significant/length(cell_outlines));
if plot_all_cells
    fprintf('Cells plotted: ALL %d cells\n', length(cells_to_plot));
else
    fprintf('Cells plotted: %d significant cells\n', length(cells_to_plot));
end
fprintf('Total place fields across plotted cells: %d\n', sum(num_place_fields));
fprintf('Mean place fields per plotted cell: %.2f\n', mean(num_place_fields));
fprintf('Median place fields per plotted cell: %.0f\n', median(num_place_fields));

valid_sic = sic_values(~isnan(sic_values));
if ~isempty(valid_sic)
    fprintf('Mean SIC (all cells): %.4f bits/spike\n', mean(valid_sic));
    sig_sic = sic_values(is_significant);
    fprintf('Mean SIC (significant cells): %.4f bits/spike\n', mean(sig_sic));
end

fprintf('\n=== Processing Complete ===\n');
if use_smoothed_outline
    fprintf('Outline style: SMOOTHED\n');
else
    fprintf('Outline style: PIXEL-PERFECT\n');
end
fprintf('Generated files:\n');
file_counter = 1;
if plot_sic_histogram
    fprintf('  %d. cells_sic_histogram_all_cells.png\n', file_counter);
    file_counter = file_counter + 1;
end
if plot_all_cells
    fprintf('  %d. cells_outline_all_regions_all.png\n', file_counter);
    file_counter = file_counter + 1;
    fprintf('  %d. cells_outline_strongest_only_all.png\n', file_counter);
    file_counter = file_counter + 1;
    fprintf('  %d. place_field_count_histogram_all.png\n', file_counter);
else
    fprintf('  %d. cells_outline_all_regions_significant.png\n', file_counter);
    file_counter = file_counter + 1;
    fprintf('  %d. cells_outline_strongest_only_significant.png\n', file_counter);
    file_counter = file_counter + 1;
    fprintf('  %d. place_field_count_histogram.png\n', file_counter);
end

fprintf('\nUsage options:\n');
fprintf('  plot_outline_overlay()                           - Default: smoothed, significant cells\n');
fprintf('  plot_outline_overlay(''pixel_perfect'')            - Pixel-perfect edges\n');
fprintf('  plot_outline_overlay(''plot_all'')                - Plot ALL cells\n');
fprintf('  plot_outline_overlay(''plot_sic'')                - Include SIC histogram\n');
fprintf('  plot_outline_overlay(''plot_all'', ''pixel_perfect'') - Combine options\n');

end