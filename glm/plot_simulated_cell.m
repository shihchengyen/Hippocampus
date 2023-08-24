function plot_simulated_cell(model)
    % Plots rate map of simulated cell from simulate_cell_FR.m
    % Specify desired parameters within simulate_cell_FR.m before calling
    % this function.
    
    % Inputs:
    % model - type of cell to be plotted. Available options:
    % 'place', 'headdirection', 'spatialview', 'ph', 'pv', 'hv', 'phv'
    
    % Code adapted from plotgridmap.m
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

    % mark out place bins under pillars
    under_pillars = nan(256, 1); j = 1;
    for i = [9:16, 25:32]
        under_pillars(j:j+15) = 40*(i-1) + [9:16, 25:32];
        j = j+16;
    end
    
    % Generate plots for each spatial variable
    if strcmp(model, 'ph') || strcmp(model, 'pv') || strcmp(model, 'hv') || strcmp(model, 'phv')
       disp('Plotting only available for single variable models.'); 
    end
    
    if strcmp(model, 'place')
        fp = figure('Name','Place plot');
        ratemap = nan(1600,1);
        for k = 1:size(ratemap,1)
            ratemap(k) = simulate_cell_FR([k NaN NaN]);
        end
        ratemap(under_pillars) = 0;

        surf(floor_x, floor_y, floor_z, flipud(reshape(ratemap(1:1600), 40, 40)'));
        alpha 1; shading flat;
        view(-35,20);
        colormap jet;
        colorbar;
    end

    if strcmp(model, 'headdirection')
        fh = figure('Name','Head direction plot');
        ratemap = nan(60,1);
        for k = 1:size(ratemap,1)
            ratemap(k) = simulate_cell_FR([NaN k NaN]);
        end

        pax = polaraxes();
        polarplot(deg2rad((0:60)*360/60), [ratemap; ratemap(1)]);
        pax.ThetaZeroLocation = 'top';
        pax.ThetaDir = 'clockwise';
    end

    if strcmp(model, 'spatialview')
        fv = figure('Name','View plot');
        ratemap = nan(5122,1);
        for k = 3:size(ratemap,1)
            ratemap(k) = simulate_cell_FR([NaN NaN k]);
        end
        ratemap(under_pillars+2) = 0;

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

end
