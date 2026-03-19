function plot_place_field_overlap()
% PLOT_PLACE_FIELD_OVERLAP - Creates a heatmap showing how many place fields overlap at each bin
% Shows density of place field coverage across the maze

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

% Extract cell outlines
fprintf('\n=== Extracting Place Field Data ===\n');
cell_outlines = {};
valid_cells = 0;

for i = 1:size(data, 1)
    rc = data{i, 1};  % region_cells
    pi = data{i, 2};  % peak_indices
    pv = data{i, 3};  % peak_values
    
    if ~isempty(rc) && iscell(rc) && ~isempty(pi) && ~isempty(pv)
        valid_cells = valid_cells + 1;
        if size(data, 2) >= 5
            cell_outlines{valid_cells} = {rc, pi, pv, data{i,4}, data{i,5}};
        elseif size(data, 2) >= 4
            cell_outlines{valid_cells} = {rc, pi, pv, data{i,4}, false};
        else
            cell_outlines{valid_cells} = {rc, pi, pv, NaN, false};
        end
    end
end

fprintf('Found %d cells with valid place field data\n', valid_cells);

if valid_cells == 0
    warning('No cells with valid data found.');
    return;
end

%% OPTIONAL: Filter to significant cells only
% Uncomment this section if you want to show only significant cells

use_significant_only = false;  % SET TO true TO FILTER

if use_significant_only
    fprintf('\n=== Filtering to Significant Cells ===\n');
    
    is_significant = false(length(cell_outlines), 1);
    
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
        
        for i = 1:length(cell_outlines)
            if i > length(cell_dirs)
                continue;
            end
            
            cell_dir = cell_dirs{i};
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
                continue;
            end
            
            % Get threshold
            if isfield(vpc_cell.data, 'critsh_sm') && ~isempty(vpc_cell.data.critsh_sm)
                threshold = prctile(vpc_cell.data.critsh_sm, 95);
                sic_value = cell_outlines{i}{4};
                if ~isnan(sic_value) && sic_value > threshold
                    is_significant(i) = true;
                end
            end
        end
    end
    
    % Filter to significant only
    significant_outlines = {};
    for i = 1:length(cell_outlines)
        if is_significant(i)
            significant_outlines{end+1} = cell_outlines{i};
        end
    end
    
    cell_outlines = significant_outlines;
    fprintf('Using %d significant cells\n', length(cell_outlines));
end

%% COUNT OVERLAPS AT EACH BIN
fprintf('\n=== Counting Place Field Overlaps ===\n');

% Initialize overlap count matrix
overlap_count = zeros(grid_size, grid_size);

% Get maze validity mask (which bins are valid, not walls)
map_2d = reshape(sample_map, grid_size, grid_size);
valid_maze = ~isnan(map_2d);

% Count overlaps for each cell's place fields
total_fields = 0;
for cell_idx = 1:length(cell_outlines)
    region_cells = cell_outlines{cell_idx}{1};
    
    if ~iscell(region_cells) || isempty(region_cells)
        continue;
    end
    
    % For each place field of this cell
    for region_idx = 1:length(region_cells)
        current_region = region_cells{region_idx};
        [rows, cols] = ind2sub([grid_size, grid_size], current_region);
        
        if length(rows) < 1
            continue;
        end
        
        % Increment count for each bin in this place field
        for i = 1:length(rows)
            r = rows(i);
            c = cols(i);
            
            % Only count if it's a valid maze location (not wall)
            if valid_maze(r, c)
                overlap_count(r, c) = overlap_count(r, c) + 1;
            end
        end
        
        total_fields = total_fields + 1;
    end
end

fprintf('Total place fields counted: %d\n', total_fields);
fprintf('Max overlap at any bin: %d place fields\n', max(overlap_count(:)));

% Set wall regions to NaN for plotting
overlap_count(~valid_maze) = NaN;

%% PLOT OVERLAP HEATMAP
fprintf('\n=== Creating Overlap Heatmap ===\n');

figure('Position', [100 100 900 800], 'Name', 'Place Field Overlap Heatmap');

% Use imagesc to plot the heatmap
h = imagesc(overlap_count);
set(h, 'AlphaData', ~isnan(overlap_count));  % Make NaN transparent

% Colormap
colormap(jet);
c = colorbar;
c.Label.String = 'Number of Overlapping Place Fields';
c.Label.FontSize = 12;
c.Label.FontWeight = 'bold';

% Axis properties
axis equal tight;
set(gca, 'YDir', 'normal');  % Flip to match standard orientation
xlabel('X Position (bins)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Y Position (bins)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Place Field Overlap Density (%d cells, %d fields)', ...
    length(cell_outlines), total_fields), 'FontSize', 14, 'FontWeight', 'bold');

% Add gridlines for pillars (optional)
hold on;
% Draw pillar outlines
pillar_coords = [8 8; 8 16; 16 16; 16 8; 8 8];  % Bottom-left pillar
plot(pillar_coords(:,1), pillar_coords(:,2), 'k-', 'LineWidth', 1.5);
pillar_coords = [24 8; 24 16; 32 16; 32 8; 24 8];  % Bottom-right pillar
plot(pillar_coords(:,1), pillar_coords(:,2), 'k-', 'LineWidth', 1.5);
pillar_coords = [8 24; 8 32; 16 32; 16 24; 8 24];  % Top-left pillar
plot(pillar_coords(:,1), pillar_coords(:,2), 'k-', 'LineWidth', 1.5);
pillar_coords = [24 24; 24 32; 32 32; 32 24; 24 24];  % Top-right pillar
plot(pillar_coords(:,1), pillar_coords(:,2), 'k-', 'LineWidth', 1.5);
hold off;

% Save figure
saveas(gcf, 'place_field_overlap_heatmap.png');
fprintf('Saved to place_field_overlap_heatmap.png\n');

%% GENERATE STATISTICS
fprintf('\n=== Overlap Statistics ===\n');

valid_counts = overlap_count(~isnan(overlap_count));
fprintf('Total maze bins: %d\n', sum(valid_maze(:)));
fprintf('Bins with at least 1 field: %d (%.1f%%)\n', sum(valid_counts > 0), ...
    100*sum(valid_counts > 0)/sum(valid_maze(:)));
fprintf('Mean overlap per bin: %.2f fields\n', mean(valid_counts));
fprintf('Median overlap per bin: %.1f fields\n', median(valid_counts));
fprintf('Max overlap at any bin: %d fields\n', max(valid_counts));
fprintf('Std dev of overlap: %.2f fields\n', std(valid_counts));

% Find hotspots (bins with highest overlap)
[max_overlap, max_idx] = max(overlap_count(:));
[max_r, max_c] = ind2sub([grid_size, grid_size], max_idx);
fprintf('\nHotspot location: bin (%d, %d) with %d overlapping fields\n', ...
    max_r, max_c, max_overlap);

%% COVERAGE HISTOGRAM
fprintf('\n=== Creating Coverage Histogram ===\n');
figure('Position', [200 200 800 500], 'Name', 'Overlap Distribution');
histogram(valid_counts, 'BinMethod', 'integers', 'FaceColor', [0.2 0.4 0.8]);
xlabel('Number of Overlapping Place Fields', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Number of Bins', 'FontSize', 12, 'FontWeight', 'bold');
title('Distribution of Place Field Overlap Across Maze', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Add statistics text
stats_text = sprintf(['Mean: %.2f\nMedian: %.1f\nMax: %d\n' ...
    'Bins covered: %d/%d (%.1f%%)'], ...
    mean(valid_counts), median(valid_counts), max(valid_counts), ...
    sum(valid_counts > 0), sum(valid_maze(:)), ...
    100*sum(valid_counts > 0)/sum(valid_maze(:)));

annotation('textbox', [0.65, 0.65, 0.25, 0.25], ...
    'String', stats_text, ...
    'FitBoxToText', 'on', ...
    'BackgroundColor', 'white', ...
    'EdgeColor', 'black', ...
    'LineWidth', 1.5, ...
    'FontSize', 11, ...
    'FontWeight', 'bold');

saveas(gcf, 'place_field_overlap_histogram.png');
fprintf('Saved to place_field_overlap_histogram.png\n');

fprintf('\n=== Processing Complete ===\n');
fprintf('Generated files:\n');
fprintf('  1. place_field_overlap_heatmap.png (2D heatmap)\n');
fprintf('  2. place_field_overlap_histogram.png (distribution)\n');

end