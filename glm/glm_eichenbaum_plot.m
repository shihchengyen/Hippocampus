function glm_eichenbaum_plot(eb_results, model, pfield_num, svfield_num)
%
%   Used to visualise the fitted results of glm_eichenbaum.
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mapping coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

place_map = flipud(reshape(1:1600, 40, 40)');

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
    case 'place'
        place_params = eb_results.params_place{pfield_num};
    case 'view'
        view_params = eb_results.params_view{svfield_num};
    case 'both'
        params = eb_results.params_joint{pfield_num, svfield_num};
        place_params = params(1:6*pfield_num);
        view_params = params(6*pfield_num + 1:end);
end

tbin_size = eb_results.tbin_size;

if strcmp(model, 'place') || strcmp(model, 'both')
    fp = figure('Name','Place plot');
    ratemap = nan(1600,1);
    for k = 1:size(ratemap,1)
        [row,col] = find(place_map == k);
        x = (floor_x(row,col) - 20) / 60;
        y = (floor_y(row,col) - 20) / 60;
        
        poiss_lambda = nan(pfield_num,1);

        for pf = 1:pfield_num
            exponent = place_params((pf-1)*6 + 1)*x + place_params((pf-1)*6 + 2)*(x^2) + ...
                place_params((pf-1)*6 + 3)*y + place_params((pf-1)*6 + 4)*(y^2) + ...
                place_params((pf-1)*6 + 5)*x*y + place_params((pf-1)*6 + 6);
            poiss_lambda(pf) = exp(exponent);
        end
        poiss_lambda = sum(poiss_lambda);
        
        ratemap(k) = poiss_lambda;
    end
    
    ratemap = ratemap / tbin_size;
    
    surf(floor_x, floor_y, floor_z, flipud(reshape(ratemap(1:1600), 40, 40)'));
    alpha 1; shading flat;
    view(-35,20);
    colormap jet;
    colorbar;
    
end
if strcmp(model, 'view') || strcmp(model, 'both')
    fv = figure('Name','View plot');
    ratemap = nan(5122,1);
    for k = 3:size(ratemap,1)
        if any(3:3+1600-1 == k) % floor
            [row,col] = find(floor == k);
            x = floor_x(row,col);
            y = floor_y(row,col);
            z = floor_z(row,col);
        elseif any(1603:1603+1600-1 == k) % ceiling
            [row,col] = find(ceiling == k);
            x = ceiling_x(row,col);
            y = ceiling_y(row,col);
            z = ceiling_z(row,col);
        elseif any(3203:3203+1280-1 == k) % walls
            [row,col] = find(walls == k);
            x = walls_x(row,col);
            y = walls_y(row,col);
            z = walls_z(row,col);
        elseif any(4483:4483+160-1 == k) % P1_BR
            [row,col] = find(P1_BR == k);
            x = P1_x(row,col);
            y = P1_y(row,col);
            z = PX_z(row,col);
        elseif any(4643:4643+160-1 == k) % P2_BL
            [row,col] = find(P2_BL == k);
            x = P2_x(row,col);
            y = P2_y(row,col);
            z = PX_z(row,col);
        elseif any(4803:4803+160-1 == k)
            [row,col] = find(P3_TR == k); % P3_TR
            x = P3_x(row,col);
            y = P3_y(row,col);
            z = PX_z(row,col);
        elseif any(4963:4963+160-1 == k) % P4_TL
            [row,col] = find(P4_TL == k);
            x = P4_x(row,col);
            y = P4_y(row,col);
            z = PX_z(row,col);
        end
        
        x = (x - 20) / 60;
        y = (y - 20) / 60;
        z = (z - 4) / 60;

        poiss_lambda = nan(svfield_num,1);
        
        for svf = 1:svfield_num
            exponent = view_params((svf-1)*10 + 1)*(x^2) + view_params((svf-1)*10 + 2)*(y^2) + view_params((svf-1)*10 + 3)*(z^2) + ...
                view_params((svf-1)*10 + 4)*x*y + view_params((svf-1)*10 + 5)*x*z + view_params((svf-1)*10 + 6)*y*z + ...
                view_params((svf-1)*10 + 7)*x + view_params((svf-1)*10 + 8)*y + view_params((svf-1)*10 + 9)*z + view_params((svf-1)*10 + 10);
            if exp(exponent) == 0
                poiss_lambda(svf) = 1e-14;
            else
                poiss_lambda(svf) = exp(exponent);
            end
        end
        poiss_lambda = sum(poiss_lambda);
        
        ratemap(k) = poiss_lambda;
    end

    ratemap = ratemap / tbin_size;

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
