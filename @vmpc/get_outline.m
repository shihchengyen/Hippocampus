function [region_cells, peak_indices, peak_values, sic_value, is_significant] = get_outline(obj, varargin)
    % Returns ALL significant connected components using half-height method
    % Peak detection: percentage of max firing rate
    % Region growing: half-height from peak to global floor
    % Display filtering: percentile threshold
    
    %USAGE 
    %plot(vpc,1,'outline') to get the default OR change the values with
    %this 
    %plot(vpc, 1, 'outline', 'outline_threshold', 0.10, 'peak_detection_pct', 0.25)
   
    Args = struct('outline_threshold', 0.75, 'Smooth', 1, 'CellIndex', 1, ...
                  'peak_detection_pct', 0.85, ...  
                  'min_component_size', 10);  % Minimum pixels for a valid component
    [Args, ~] = getOptArgs(varargin, Args);

    if Args.Smooth
        map_1d = obj.data.maps_adsm(Args.CellIndex, :);
    else
        map_1d = obj.data.maps_raw(Args.CellIndex, :);
    end

    grid_size = sqrt(length(map_1d));
    map_2d = reshape(map_1d, grid_size, grid_size);

    % Remove NaN values for statistics
    valid_values = map_2d(~isnan(map_2d));
    
    if isempty(valid_values)
        region_cells = {};
        peak_indices = [];
        peak_values = [];
        sic_value = NaN;
        is_significant = false;
        return;
    end
    
    % Calculate statistics
    map_mean = mean(valid_values);
    map_std = std(valid_values);
    map_max = max(valid_values);
    map_min = min(valid_values);  % Global minimum (floor)
    
    % Peak detection threshold - percentage of max firing rate
    peak_detection_threshold = Args.peak_detection_pct * map_max;
    
    % Display threshold (for filtering which regions to show)
    display_threshold = prctile(valid_values, Args.outline_threshold * 100);
    
    fprintf('Cell %d - Map stats: mean=%.2f, std=%.2f, max=%.2f, min=%.2f\n', ...
        Args.CellIndex, map_mean, map_std, map_max, map_min);
    fprintf('  Peak Detection (%.0f%% of max firing rate): %.2f Hz\n', ...
        Args.peak_detection_pct * 100, peak_detection_threshold);
    fprintf('  Display Filter (%.0fth percentile): %.2f Hz\n', ...
        Args.outline_threshold * 100, display_threshold);
    fprintf('  Half-height baseline (global floor): %.2f Hz\n', map_min);
    
    % Step 1: Find ALL local maxima above peak detection threshold
    all_peaks = [];
    
    % Find all pixels above peak detection threshold
    above_threshold = (map_2d >= peak_detection_threshold) & ~isnan(map_2d);
    
    % For each pixel above threshold, check if it's a local maximum
    for r = 1:grid_size
        for c = 1:grid_size
            if ~above_threshold(r, c)
                continue;
            end
            
            current_val = map_2d(r, c);
            
            % Check 8-connected neighbors
            is_peak = true;
            for dr = -1:1
                for dc = -1:1
                    if dr == 0 && dc == 0
                        continue;
                    end
                    
                    nr = r + dr;
                    nc = c + dc;
                    
                    if nr >= 1 && nr <= grid_size && nc >= 1 && nc <= grid_size
                        if ~isnan(map_2d(nr, nc)) && map_2d(nr, nc) > current_val
                            is_peak = false;
                            break;
                        end
                    end
                end
                if ~is_peak
                    break;
                end
            end
            
            if is_peak
                peak_idx = sub2ind([grid_size, grid_size], r, c);
                all_peaks = [all_peaks; r, c, current_val, peak_idx];
            end
        end
    end
    
    num_peaks = size(all_peaks, 1);
    fprintf('  Found %d peaks above %.0f%% of max firing rate\n', num_peaks, Args.peak_detection_pct * 100);
    
    % Get SIC value and significance for this cell
    if isfield(obj.data, 'SIC_adsm') && Args.CellIndex <= size(obj.data.SIC_adsm, 1)
        sic_value = obj.data.SIC_adsm(Args.CellIndex, 1);
        
        % Check if SIC is significant using critsh_sm threshold
        if isfield(obj.data, 'critsh_sm')
            % Get 95th percentile of critsh_sm (global threshold)
            shuffle_threshold = prctile(obj.data.critsh_sm, 95);
            is_significant = (sic_value > shuffle_threshold);
            
            fprintf('  SIC: %.4f bits/spike\n', sic_value);
            fprintf('  Shuffle threshold (95%%): %.4f\n', shuffle_threshold);
            if is_significant
                fprintf('  *** SIGNIFICANT place cell (p < 0.05) ***\n');
            else
                fprintf('  Not significant (p > 0.05)\n');
            end
        else
            is_significant = false;
            fprintf('  SIC: %.4f bits/spike (no shuffle data)\n', sic_value);
        end
    else
        sic_value = NaN;
        is_significant = false;
        fprintf('  SIC: not available\n');
    end
    
    if num_peaks == 0
        region_cells = {};
        peak_indices = [];
        peak_values = [];
        return;
    end
    
    % Step 2: For each peak, grow region using HALF-HEIGHT method
    % Use 4-connectivity and simple flood fill from each peak
    all_components = {};
    component_peaks = {};
    
    % Track globally assigned pixels to prevent overlap
    globally_assigned = false(grid_size, grid_size);
    
    % Sort peaks by height (process strongest first)
    [~, peak_order] = sort(all_peaks(:, 3), 'descend');
    
    for p_idx = 1:num_peaks
        p = peak_order(p_idx);
        peak_row = all_peaks(p, 1);
        peak_col = all_peaks(p, 2);
        peak_val = all_peaks(p, 3);
        
        % Skip if this peak location is already claimed
        if globally_assigned(peak_row, peak_col)
            fprintf('    Peak %d at (%.0f,%.0f) already claimed, skipping\n', ...
                p, peak_row, peak_col);
            continue;
        end
        
        % Calculate half-height threshold for THIS peak
        half_height_threshold = peak_val - (peak_val - map_min) / 2;
        
        fprintf('    Peak %d at (%.0f,%.0f) = %.2f Hz, half-height: %.2f Hz\n', ...
            p, peak_row, peak_col, peak_val, half_height_threshold);
        
        % Flood fill from this peak using 4-connectivity (more restrictive)
        region_mask = false(grid_size, grid_size);
        region_mask(peak_row, peak_col) = true;
        
        queue = [peak_row, peak_col];
        head = 1;
        
        while head <= size(queue, 1)
            r = queue(head, 1);
            c = queue(head, 2);
            head = head + 1;
            
            % Only 4 neighbors (up, down, left, right) - more restrictive than 8
            neighbors = [r-1, c; r+1, c; r, c-1; r, c+1];
            
            for i = 1:4
                nr = neighbors(i, 1);
                nc = neighbors(i, 2);
                
                if nr >= 1 && nr <= grid_size && nc >= 1 && nc <= grid_size
                    % Check: not visited, above half-height, not NaN, not globally assigned
                    if ~region_mask(nr, nc) && ...
                       ~globally_assigned(nr, nc) && ...
                       ~isnan(map_2d(nr, nc)) && ...
                       map_2d(nr, nc) >= half_height_threshold
                        
                        region_mask(nr, nc) = true;
                        queue = [queue; nr, nc];
                    end
                end
            end
        end
        
        region_size = sum(region_mask(:));
        
        if region_size >= Args.min_component_size
            % Mark these pixels as globally assigned
            globally_assigned = globally_assigned | region_mask;
            
            all_components{end+1} = region_mask;
            component_peaks{end+1} = all_peaks(p, :);
            
            fprintf('      Component %d: %d pixels\n', length(all_components), region_size);
        else
            fprintf('      Skipped: region too small (%d pixels)\n', region_size);
        end
    end
    
    % Step 3: Return ALL components, but FILTER by display threshold
    if isempty(all_components)
        region_cells = {};
        peak_indices = [];
        peak_values = [];
        fprintf('  No components found\n\n');
        return;
    end
    
    % Sort components by peak value (strongest first)
    component_peak_vals = zeros(length(all_components), 1);
    for i = 1:length(all_components)
        component_peak_vals(i) = component_peaks{i}(3);  % Peak value
    end
    [~, sort_idx] = sort(component_peak_vals, 'descend');
    
    fprintf('  Found %d connected components\n', length(all_components));
    
    % Return components that pass display threshold
    region_cells = {};
    peak_indices = [];
    peak_values = [];
    num_filtered = 0;
    
    for i = 1:length(all_components)
        idx = sort_idx(i);
        component_mask = all_components{idx};
        region_idx = find(component_mask);
        
        % Get the peak info for this component
        primary_peak = component_peaks{idx};
        peak_value = primary_peak(3);
        
        % Check if peak exceeds display threshold
        if peak_value >= display_threshold
            % Store this region
            region_cells{end+1} = region_idx;
            peak_indices = [peak_indices; primary_peak(4)];
            peak_values = [peak_values; peak_value];
            
            region_values = map_2d(region_idx);
            fprintf('  Component %d: %d pixels (%.2f-%.2f Hz), peak at (%.0f,%.0f) = %.2f Hz [DISPLAYED]\n', ...
                i, length(region_idx), min(region_values), max(region_values), ...
                primary_peak(1), primary_peak(2), peak_value);
        else
            num_filtered = num_filtered + 1;
            fprintf('  Component %d: peak = %.2f Hz < %.2f Hz [FILTERED OUT]\n', ...
                i, peak_value, display_threshold);
        end
    end
    
    fprintf('  Total: %d regions displayed (peak > %.0fth percentile = %.2f Hz)\n', ...
        length(region_cells), Args.outline_threshold * 100, display_threshold);
    fprintf('         %d regions filtered out (peak below %.0fth percentile)\n\n', ...
        num_filtered, Args.outline_threshold * 100);
end