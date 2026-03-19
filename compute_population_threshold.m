function population_threshold = compute_population_threshold()
% COMPUTE_POPULATION_THRESHOLD Computes population-level SIC threshold
% and identifies significant cells
% 
% This function:
% 1. Loads all cells from place_selective_cells.txt
% 2. Extracts shuffle distributions from each cell's vmpc.mat
% 3. Pools all shuffle values together
% 4. Computes the 95th percentile as the population threshold
% 5. Loads actual SIC values and identifies significant cells
% 6. Creates histogram showing significant vs non-significant cells
%
% Output:
%   population_threshold - Single threshold value for the entire population

fprintf('\n========================================\n');
fprintf('COMPUTING POPULATION THRESHOLD\n');
fprintf('========================================\n\n');

% Check if place_selective_cells.txt exists
if ~exist('place_selective_cells.txt', 'file')
    error('place_selective_cells.txt not found in current directory');
end

% Load cell directories from place_selective_cells.txt
fid = fopen('place_selective_cells.txt', 'r');
cell_dirs = {};
while ~feof(fid)
    line = fgetl(fid);
    if ischar(line) && ~isempty(line)
        cell_dirs{end+1} = line;
    end
end
fclose(fid);

fprintf('Found %d cell directories in place_selective_cells.txt\n\n', length(cell_dirs));

if isempty(cell_dirs)
    error('No cell directories found in place_selective_cells.txt');
end

% Initialize arrays
all_shuffled_sics = [];
actual_sics = zeros(length(cell_dirs), 1);
successful_cells = 0;
failed_cells = 0;

% Loop through each cell and extract shuffle distributions AND actual SIC
fprintf('Loading data from each cell...\n');
fprintf('-------------------------------------------\n');

for i = 1:length(cell_dirs)
    cell_dir = cell_dirs{i};
    
    % Try multiple possible locations for vmpc.mat
    possible_paths = {
        fullfile(cell_dir, 'FiltAll', '1px', 'vmpc.mat'),
        fullfile(cell_dir, 'FiltAll', 'vmpc.mat'),
        fullfile(cell_dir, 'vmpc.mat')
    };
    
    vpc_cell = [];
    loaded_path = '';
    
    % Try to load from each possible path
    for p = 1:length(possible_paths)
        if exist(possible_paths{p}, 'file')
            try
                vpc_temp = load(possible_paths{p});
                if isfield(vpc_temp, 'vmp')
                    vpc_cell = vpc_temp.vmp;
                    loaded_path = possible_paths{p};
                    break;
                elseif isfield(vpc_temp, 'vpc')
                    vpc_cell = vpc_temp.vpc;
                    loaded_path = possible_paths{p};
                    break;
                end
            catch
                continue;
            end
        end
    end
    
    % Extract shuffle distribution AND actual SIC if vmpc was loaded successfully
    if ~isempty(vpc_cell)
        shuffle_values = [];
        
        % Get shuffle distribution
        if isfield(vpc_cell.data, 'critsh_sm') && ~isempty(vpc_cell.data.critsh_sm)
            shuffle_values = vpc_cell.data.critsh_sm(:);
        elseif isfield(vpc_cell.data, 'SICsh') && ~isempty(vpc_cell.data.SICsh)
            shuffle_values = vpc_cell.data.SICsh(:);
        end
        
        % Get actual SIC value
        if isfield(vpc_cell.data, 'SIC') && ~isempty(vpc_cell.data.SIC)
            actual_sics(i) = vpc_cell.data.SIC;
        elseif isfield(vpc_cell.data, 'sic') && ~isempty(vpc_cell.data.sic)
            actual_sics(i) = vpc_cell.data.sic;
        end
        
        if ~isempty(shuffle_values)
            all_shuffled_sics = [all_shuffled_sics; shuffle_values];
            successful_cells = successful_cells + 1;
            fprintf('Cell %3d: ? SIC=%.4f, %d shuffle values loaded\n', ...
                i, actual_sics(i), length(shuffle_values));
        else
            failed_cells = failed_cells + 1;
            fprintf('Cell %3d: ? No shuffle data found\n', i);
        end
    else
        failed_cells = failed_cells + 1;
        fprintf('Cell %3d: ? Could not load vmpc.mat\n', i);
        actual_sics(i) = NaN;
    end
end

fprintf('\n-------------------------------------------\n');
fprintf('Loading Summary:\n');
fprintf('  Successfully loaded: %d cells\n', successful_cells);
fprintf('  Failed to load: %d cells\n', failed_cells);
fprintf('  Total shuffle values: %d\n', length(all_shuffled_sics));
fprintf('-------------------------------------------\n\n');

if isempty(all_shuffled_sics)
    error('No shuffle data found. Cannot compute population threshold.');
end

% Compute the population threshold (95th percentile)
population_threshold = prctile(all_shuffled_sics, 95);

% Identify significant cells
valid_sics = ~isnan(actual_sics);
is_significant = (actual_sics > population_threshold) & valid_sics;
num_significant = sum(is_significant);

fprintf('========================================\n');
fprintf('POPULATION THRESHOLD RESULTS\n');
fprintf('========================================\n');
fprintf('Population threshold (95th percentile): %.6f bits/spike\n', population_threshold);
fprintf('Based on %d shuffle values from %d cells\n\n', length(all_shuffled_sics), successful_cells);
fprintf('SIGNIFICANCE RESULTS:\n');
fprintf('  Total cells: %d\n', sum(valid_sics));
fprintf('  Significant cells: %d (%.1f%%)\n', num_significant, 100*num_significant/sum(valid_sics));
fprintf('  Non-significant cells: %d (%.1f%%)\n', sum(valid_sics)-num_significant, 100*(sum(valid_sics)-num_significant)/sum(valid_sics));
fprintf('========================================\n\n');

% Create histogram of ACTUAL SIC values with threshold line
figure('Position', [100 100 1000 600], 'Name', 'Cell Significance Based on Population Threshold');

% Prepare data
cell_numbers = 1:length(actual_sics);
valid_idx = valid_sics;
sic_vals = actual_sics(valid_idx);
sig_idx = is_significant(valid_idx);

% Plot non-significant cells in gray
non_sig_cells = cell_numbers(valid_idx);
non_sig_sics = sic_vals;
non_sig_cells(sig_idx) = [];
non_sig_sics(sig_idx) = [];

if ~isempty(non_sig_cells)
    bar(non_sig_cells, non_sig_sics, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'k');
end

hold on;

% Plot significant cells in green
sig_cells = cell_numbers(valid_idx);
sig_sics = sic_vals;
sig_cells = sig_cells(sig_idx);
sig_sics = sig_sics(sig_idx);

if ~isempty(sig_cells)
    bar(sig_cells, sig_sics, 'FaceColor', [0.2 0.8 0.2], 'EdgeColor', 'k');
end

hold on;

% Add population threshold line
x_lim = xlim;
plot([x_lim(1) x_lim(2)], [population_threshold population_threshold], ...
    'r-', 'LineWidth', 3);
plot([x_lim(1) x_lim(2)], [population_threshold population_threshold], ...
    'r--', 'LineWidth', 1.5);

% Add asterisks above significant cells
if any(sig_idx)
    sig_cell_nums = cell_numbers(valid_idx);
    sig_cell_nums = sig_cell_nums(sig_idx);
    sig_sic_vals = sic_vals(sig_idx);
    plot(sig_cell_nums, sig_sic_vals + 0.01, 'g*', 'MarkerSize', 12, 'LineWidth', 2);
end

hold off;

xlabel('Cell Number', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('SIC (bits/spike)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Cell SIC Values | Red Line = Population Threshold (%.4f)', population_threshold), ...
    'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Add statistics box
stats_text = sprintf(['Total: %d cells\nSignificant: %d (%.1f%%)\n' ...
    'Non-significant: %d (%.1f%%)\nThreshold: %.4f'], ...
    sum(valid_sics), num_significant, 100*num_significant/sum(valid_sics), ...
    sum(valid_sics)-num_significant, 100*(sum(valid_sics)-num_significant)/sum(valid_sics), ...
    population_threshold);

annotation('textbox', [0.15, 0.70, 0.25, 0.2], ...
    'String', stats_text, ...
    'FitBoxToText', 'on', ...
    'BackgroundColor', 'white', ...
    'EdgeColor', 'black', ...
    'LineWidth', 1.5, ...
    'FontSize', 10, ...
    'FontWeight', 'bold');

% Save figure
saveas(gcf, 'population_threshold_significance.png');
fprintf('Saved figure to: population_threshold_significance.png\n\n');

% Save results to file
fid = fopen('population_threshold_results.txt', 'w');
fprintf(fid, 'POPULATION THRESHOLD ANALYSIS\n');
fprintf(fid, '================================\n\n');
fprintf(fid, 'Population Threshold (95th percentile): %.6f bits/spike\n', population_threshold);
fprintf(fid, 'Based on %d shuffle values from %d cells\n\n', length(all_shuffled_sics), successful_cells);
fprintf(fid, 'SIGNIFICANCE RESULTS:\n');
fprintf(fid, 'Total cells: %d\n', sum(valid_sics));
fprintf(fid, 'Significant cells: %d (%.1f%%)\n', num_significant, 100*num_significant/sum(valid_sics));
fprintf(fid, 'Non-significant cells: %d (%.1f%%)\n\n', sum(valid_sics)-num_significant, 100*(sum(valid_sics)-num_significant)/sum(valid_sics));
fprintf(fid, 'SIGNIFICANT CELLS:\n');
fprintf(fid, 'Cell#\tSIC\n');
for i = 1:length(actual_sics)
    if is_significant(i)
        fprintf(fid, '%d\t%.6f\n', i, actual_sics(i));
    end
end
fclose(fid);

fprintf('Saved results to: population_threshold_results.txt\n');
fprintf('\n========================================\n');
fprintf('DONE!\n');
fprintf('========================================\n');

end