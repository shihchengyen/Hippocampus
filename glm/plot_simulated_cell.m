function plot_simulated_cell()
    % Plots rate map of simulated cell from simulate_cell_FR.m
    % Specify desired parameters within simulate_cell_FR.m before calling
    % this function.
    
    % read type of cell simulated by simulate_cell_FR.m
    [~, cell_params] = simulate_cell_FR([NaN, NaN, NaN]);
    model = cell_params.cell_type;
    
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

    % mark out place bins under pillars
    under_pillars = nan(4*pillar_width^2, 1); j = 1;
    for i = [pillar_width+1:2*pillar_width, 3*pillar_width+1:4*pillar_width]
        under_pillars(j:j+2*pillar_width-1) = floor_width*(i-1) + [pillar_width+1:2*pillar_width, 3*pillar_width+1:4*pillar_width];
        j = j+2*pillar_width;
    end
    
    % Generate plots for each spatial variable
    if strcmp(model, 'place') || strcmp(model, 'ph') || strcmp(model, 'pv') || strcmp(model, 'phv')
        fp = figure('Name','Place plot');
        place_bins = [(1:num_place_bins)', nan(num_place_bins,1), nan(num_place_bins,1)]; 
        ratemap = simulate_cell_FR(place_bins);
        ratemap(under_pillars) = 0;

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
    end

    if strcmp(model, 'headdirection') || strcmp(model, 'ph') || strcmp(model, 'hv') || strcmp(model, 'phv')
        fh = figure('Name','Head direction plot');
        hd_bins = [nan(num_hd_bins,1), (1:num_hd_bins)', nan(num_hd_bins,1)]; 
        ratemap = simulate_cell_FR(hd_bins);

        pax = polaraxes();
        polarplot(deg2rad((0:num_hd_bins)*360/num_hd_bins), [ratemap; ratemap(1)]);
        pax.ThetaZeroLocation = 'top';
        pax.ThetaDir = 'clockwise';
    end

    if strcmp(model, 'spatialview') || strcmp(model, 'pv') || strcmp(model, 'hv') || strcmp(model, 'phv')
        fv = figure('Name','View plot');
        view_bins = [nan(num_view_bins,1), nan(num_view_bins,1), (1:num_view_bins)']; 
        ratemap = simulate_cell_FR(view_bins);
        ratemap(under_pillars+viewbin_offset) = 0;

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
    end

end
