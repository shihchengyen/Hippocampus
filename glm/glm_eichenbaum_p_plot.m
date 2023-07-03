function glm_eichenbaum_p_plot(eb_results, model, pfield_num, svfield_num)
%
%   Used to visualise the fitted results of glm_eichenbaum.
%

padding = eb_results.padding;
tbin_size = eb_results.tbin_size;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mapping coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Code adapted from plotgridmap.m

floor_x = repmat(0:40, 41, 1);
floor_y = flipud(repmat([0:40]', 1, 41));
floor_z = zeros(41,41);

sv_floor_x = floor_x;
sv_floor_y = floor_y;

ceiling_x = floor_x;
ceiling_y = floor_y;
ceiling_z = 8.*ones(41,41);

sv_ceiling_x = sv_floor_x;
sv_ceiling_y = flipud(flipud(repmat([48:88]', 1, 41)));

walls_x = repmat([0.*ones(1,40) 0:39 40.*ones(1,40) 40:-1:0], 9, 1); % left wall, top wall, right wall, bottom wall. Order is clockwise, top to bottom.
walls_y = repmat([0:39 40.*ones(1,40) 40:-1:1 0.*ones(1,41)], 9, 1);
walls_z = repmat([8:-1:0]', 1, 40*4 + 1);

sv_walls_x = repmat([-40:-1 0:39 40:79 80:120], 9, 1);
sv_walls_y = repmat([48:-1:40]', 1, 40*4 + 1);

P1_x = repmat([24.*ones(1,8) 24:31 32.*ones(1,8) 32:-1:24], 6, 1); % clockwise starting from left face
P1_y = repmat([8:15 16.*ones(1,8) 16:-1:9 8.*ones(1,9)], 6, 1);
PX_z = repmat([5:-1:0]', 1, 8*4 + 1);

sv_P1_x = repmat([-120:-88], 6, 1);
sv_P1_y = repmat([5:-1:0]', 1, 8*4 + 1);

P2_x = repmat([8.*ones(1,8) 8:15 16.*ones(1,8) 16:-1:8], 6, 1); % clockwise starting from left face
P2_y = P1_y;

sv_P2_x = sv_P1_x;
sv_P2_y = sv_P1_y + 20;

P3_x = P1_x; % clockwise starting from left face
P3_y = repmat([24:31 32.*ones(1,8) 32:-1:25 24.*ones(1,9)], 6, 1);

sv_P3_x = sv_P1_x;
sv_P3_y = sv_P1_y + 40;

P4_x = P2_x; % clockwise starting from left face
P4_y = P3_y;

sv_P4_x = sv_P1_x;
sv_P4_y = sv_P1_y + 60;

% Links: if the data has coordinates that belong in a Link, it will be
% duplicated and matched to the corresponding coordinates with the same
% indices in its Pad.
% For links: always arrange so that x is increasing, y is decreasing as
% row/col indices increase

link_ceilingToWallL_x = sv_ceiling_x(1:end-1, 1:padding);
link_ceilingToWallL_y = sv_ceiling_y(1:end-1, 1:padding);
pad_ceilingToWallL_x = rot90(repmat([-40:0], padding+1, 1), -1);
pad_ceilingToWallL_y = rot90(repmat([48+padding:-1:48]', 1, 41), -1);

link_ceilingToWallR_x = sv_ceiling_x(1:end-1, end-padding+1:end);
link_ceilingToWallR_y = sv_ceiling_y(1:end-1, end-padding+1:end);
pad_ceilingToWallR_x = rot90(repmat([40:80], padding+1, 1), +1);
pad_ceilingToWallR_y = rot90(repmat([48+padding:-1:48]', 1, 41), +1);

link_ceilingToWallB_x = sv_ceiling_x(1:padding, 1:end-1);
link_ceilingToWallB_y = sv_ceiling_y(1:padding, 1:end-1);
pad_ceilingToWallB_x = rot90(repmat([80:120], padding+1, 1), +2);
pad_ceilingToWallB_y = rot90(repmat([48+padding:-1:48]', 1, 41), +2);

link_wallBToWallL_x = sv_walls_x(1:end-1, end-padding+1:end);
link_wallBToWallL_y = sv_walls_y(1:end-1, end-padding+1:end);
pad_wallBToWallL_x = repmat([-40-padding:-40], 9, 1);
pad_wallBToWallL_y = repmat([48:-1:40]', 1, padding+1);

% % % % % % % %

link_floorToWallL_x = sv_floor_x(1:end-1, 1:padding);
link_floorToWallL_y = sv_floor_y(1:end-1, 1:padding);
pad_floorToWallL_x = rot90(repmat([-40:0], padding+1, 1), +1);
pad_floorToWallL_y = rot90(repmat([40:-1:40-padding]', 1, 41), +1);

link_floorToWallR_x = sv_floor_x(1:end-1, end-padding+1:end);
link_floorToWallR_y = sv_floor_y(1:end-1, end-padding+1:end);
pad_floorToWallR_x = rot90(repmat([40:80], padding+1, 1), -1);
pad_floorToWallR_y = rot90(repmat([40:-1:40-padding]', 1, 41), -1);

link_floorToWallB_x = sv_floor_x(end-padding+1:end, 1:end-1);
link_floorToWallB_y = sv_floor_y(end-padding+1:end, 1:end-1);
pad_floorToWallB_x = rot90(repmat([80:120], padding+1, 1), +2);
pad_floorToWallB_y = rot90(repmat([40:-1:40-padding]', 1, 41), +2);

% % % % % % % %

link_floorToP1_x = nan(padding+8+padding, padding+8+padding);
link_floorToP1_x(:, padding+1:end-padding) = repmat([24:31], padding+8+padding, 1); % vertical strip
link_floorToP1_x(padding+1:end-padding, :) = repmat([24-padding:31+padding], 8, 1); % horizontal strip
link_floorToP1_x(padding+1:end-padding, padding+1:end-padding) = nan; % clear center
link_floorToP1_y = nan(padding+8+padding, padding+8+padding);
link_floorToP1_y(:, padding+1:end-padding) = repmat([15+padding:-1:8-padding]', 1, 8); % vertical strip
link_floorToP1_y(padding+1:end-padding, :) = repmat([15:-1:8]', 1, padding+8+padding); % horizontal strip
link_floorToP1_y(padding+1:end-padding, padding+1:end-padding) = nan; % clear center

pad_floorToP1_x = nan(padding+8+padding+1, padding+8+padding+1);
pad_floorToP1_x(padding+1:end-padding-1, 1:padding) = rot90(repmat([-120:-113], padding, 1), -1); % left face
pad_floorToP1_x(1:padding, padding+1:end-padding-1) = rot90(repmat([-112:-105], padding, 1), +2); % top face
pad_floorToP1_x(padding+1:end-padding-1, end-1-padding+1:end) = rot90(repmat([-104:-97], padding+1, 1), +1); % right face
pad_floorToP1_x(end-1-padding+1:end, padding+1:end-padding-1) = rot90(repmat([-96:-89], padding+1, 1), 0); % bottom face
pad_floorToP1_y = nan(padding+8+padding+1, padding+8+padding+1);
pad_floorToP1_y(padding+1:end-padding-1, 1:padding) = rot90(repmat([sv_P1_y(end,1):-1:sv_P1_y(end,1)-padding+1]', 1, 8), -1); % left face
pad_floorToP1_y(1:padding, padding+1:end-padding-1) = rot90(repmat([sv_P1_y(end,1):-1:sv_P1_y(end,1)-padding+1]', 1, 8), +2); % top face
pad_floorToP1_y(padding+1:end-padding-1, end-1-padding+1:end) = rot90(repmat([sv_P1_y(end,1):-1:sv_P1_y(end,1)-padding]', 1, 8), +1); % right face
pad_floorToP1_y(end-1-padding+1:end, padding+1:end-padding-1) = rot90(repmat([sv_P1_y(end,1):-1:sv_P1_y(end,1)-padding]', 1, 8), 0); % bottom face

% % % % % % % %

link_floorToP2_x = nan(padding+8+padding, padding+8+padding);
link_floorToP2_x(:, padding+1:end-padding) = repmat([8:15], padding+8+padding, 1); % vertical strip
link_floorToP2_x(padding+1:end-padding, :) = repmat([8-padding:15+padding], 8, 1); % horizontal strip
link_floorToP2_x(padding+1:end-padding, padding+1:end-padding) = nan; % clear center
link_floorToP2_y = link_floorToP1_y;

pad_floorToP2_x = pad_floorToP1_x;
pad_floorToP2_y = nan(padding+8+padding+1, padding+8+padding+1);
pad_floorToP2_y(padding+1:end-padding-1, 1:padding) = rot90(repmat([sv_P2_y(end,1):-1:sv_P2_y(end,1)-padding+1]', 1, 8), -1); % left face
pad_floorToP2_y(1:padding, padding+1:end-padding-1) = rot90(repmat([sv_P2_y(end,1):-1:sv_P2_y(end,1)-padding+1]', 1, 8), +2); % top face
pad_floorToP2_y(padding+1:end-padding-1, end-1-padding+1:end) = rot90(repmat([sv_P2_y(end,1):-1:sv_P2_y(end,1)-padding]', 1, 8), +1); % right face
pad_floorToP2_y(end-1-padding+1:end, padding+1:end-padding-1) = rot90(repmat([sv_P2_y(end,1):-1:sv_P2_y(end,1)-padding]', 1, 8), 0); % bottom face

% % % % % % % %

link_floorToP3_x = link_floorToP1_x;
link_floorToP3_y = nan(padding+8+padding, padding+8+padding);
link_floorToP3_y(:, padding+1:end-padding) = repmat([31+padding:-1:24-padding]', 1, 8); % vertical strip
link_floorToP3_y(padding+1:end-padding, :) = repmat([31:-1:24]', 1, padding+8+padding); % horizontal strip
link_floorToP3_y(padding+1:end-padding, padding+1:end-padding) = nan; % clear center

pad_floorToP3_x = pad_floorToP1_x;
pad_floorToP3_y = nan(padding+8+padding+1, padding+8+padding+1);
pad_floorToP3_y(padding+1:end-padding-1, 1:padding) = rot90(repmat([sv_P3_y(end,1):-1:sv_P3_y(end,1)-padding+1]', 1, 8), -1); % left face
pad_floorToP3_y(1:padding, padding+1:end-padding-1) = rot90(repmat([sv_P3_y(end,1):-1:sv_P3_y(end,1)-padding+1]', 1, 8), +2); % top face
pad_floorToP3_y(padding+1:end-padding-1, end-1-padding+1:end) = rot90(repmat([sv_P3_y(end,1):-1:sv_P3_y(end,1)-padding]', 1, 8), +1); % right face
pad_floorToP3_y(end-1-padding+1:end, padding+1:end-padding-1) = rot90(repmat([sv_P3_y(end,1):-1:sv_P3_y(end,1)-padding]', 1, 8), 0); % bottom face

% % % % % % % %

link_floorToP4_x = link_floorToP2_x;
link_floorToP4_y = link_floorToP3_y;

pad_floorToP4_x = pad_floorToP1_x;
pad_floorToP4_y = nan(padding+8+padding+1, padding+8+padding+1);
pad_floorToP4_y(padding+1:end-padding-1, 1:padding) = rot90(repmat([sv_P4_y(end,1):-1:sv_P4_y(end,1)-padding+1]', 1, 8), -1); % left face
pad_floorToP4_y(1:padding, padding+1:end-padding-1) = rot90(repmat([sv_P4_y(end,1):-1:sv_P4_y(end,1)-padding+1]', 1, 8), +2); % top face
pad_floorToP4_y(padding+1:end-padding-1, end-1-padding+1:end) = rot90(repmat([sv_P4_y(end,1):-1:sv_P4_y(end,1)-padding]', 1, 8), +1); % right face
pad_floorToP4_y(end-1-padding+1:end, padding+1:end-padding-1) = rot90(repmat([sv_P4_y(end,1):-1:sv_P4_y(end,1)-padding]', 1, 8), 0); % bottom face

% % % % % % % %

link_P1BToP1L_x = sv_P1_x(1:end-1, end-padding:end-1);
link_P1BToP1L_y = sv_P1_y(1:end-1, end-padding:end-1);
pad_P1BToP1L_x = repmat([-120-padding:-120], 6, 1);
pad_P1BToP1L_y = sv_P1_y(1:end, end-padding:end);

link_P2BToP2L_x = sv_P2_x(1:end-1, end-padding:end-1);
link_P2BToP2L_y = sv_P2_y(1:end-1, end-padding:end-1);
pad_P2BToP2L_x = pad_P1BToP1L_x;
pad_P2BToP2L_y = sv_P2_y(1:end, end-padding:end);

link_P3BToP3L_x = sv_P3_x(1:end-1, end-padding:end-1);
link_P3BToP3L_y = sv_P3_y(1:end-1, end-padding:end-1);
pad_P3BToP3L_x = pad_P1BToP1L_x;
pad_P3BToP3L_y = sv_P3_y(1:end, end-padding:end);

link_P4BToP4L_x = sv_P4_x(1:end-1, end-padding:end-1);
link_P4BToP4L_y = sv_P4_y(1:end-1, end-padding:end-1);
pad_P4BToP4L_x = pad_P1BToP1L_x;
pad_P4BToP4L_y = sv_P4_y(1:end, end-padding:end);

% % % % % % % %

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



if strcmp(model, 'place') || strcmp(model, 'both')
    fp = figure('Name','Place plot');
    ratemap = nan(size(floor_x)-1);
    for row = 1:size(ratemap,1)
        for col = 1:size(ratemap,2)
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

            ratemap(row,col) = poiss_lambda;
        end
    end
    
    ratemap = ratemap / tbin_size;
    
    surf(floor_x, floor_y, floor_z, ratemap);
    alpha 1; shading flat;
    view(-35,20);
    colormap jet;
    colorbar;
    
end
if strcmp(model, 'view') || strcmp(model, 'both')
    fv = figure('Name','View plot');
    
    ratemap = cell(0);
    
    ratemap{1} = nan(size(sv_floor_x)-1); coordmap_x{1} = sv_floor_x; coordmap_y{1} = sv_floor_y;
    ratemap{2} = nan(size(sv_ceiling_x)-1); coordmap_x{2} = sv_ceiling_x; coordmap_y{2} = sv_ceiling_y;
    ratemap{3} = nan(size(sv_walls_x)-1); coordmap_x{3} = sv_walls_x; coordmap_y{3} = sv_walls_y;
    ratemap{4} = nan(size(sv_P1_x)-1); coordmap_x{4} = sv_P1_x; coordmap_y{4} = sv_P1_y;
    ratemap{5} = nan(size(sv_P2_x)-1); coordmap_x{5} = sv_P2_x; coordmap_y{5} = sv_P2_y;
    ratemap{6} = nan(size(sv_P3_x)-1); coordmap_x{6} = sv_P3_x; coordmap_y{6} = sv_P3_y;
    ratemap{7} = nan(size(sv_P4_x)-1); coordmap_x{7} = sv_P4_x; coordmap_y{7} = sv_P4_y;
    
    ratemap{8} = nan(size(pad_ceilingToWallL_x)-1); coordmap_x{8} = pad_ceilingToWallL_x; coordmap_y{8} = pad_ceilingToWallL_y;
    ratemap{9} = nan(size(pad_ceilingToWallR_x)-1); coordmap_x{9} = pad_ceilingToWallR_x; coordmap_y{9} = pad_ceilingToWallR_y;
    ratemap{10} = nan(size(pad_ceilingToWallB_x)-1); coordmap_x{10} = pad_ceilingToWallB_x; coordmap_y{10} = pad_ceilingToWallB_y;
    
    ratemap{11} = nan(size(pad_wallBToWallL_x)-1); coordmap_x{11} = pad_wallBToWallL_x; coordmap_y{11} = pad_wallBToWallL_y;
    
    ratemap{12} = nan(size(pad_floorToWallL_x)-1); coordmap_x{12} = pad_floorToWallL_x; coordmap_y{12} = pad_floorToWallL_y;
    ratemap{13} = nan(size(pad_floorToWallR_x)-1); coordmap_x{13} = pad_floorToWallR_x; coordmap_y{13} = pad_floorToWallR_y;
    ratemap{14} = nan(size(pad_floorToWallB_x)-1); coordmap_x{14} = pad_floorToWallB_x; coordmap_y{14} = pad_floorToWallB_y;
    
    ratemap{15} = nan(size(pad_floorToP1_x)-1); coordmap_x{15} = pad_floorToP1_x; coordmap_y{15} = pad_floorToP1_y;
    ratemap{16} = nan(size(pad_floorToP2_x)-1); coordmap_x{16} = pad_floorToP2_x; coordmap_y{16} = pad_floorToP2_y;
    ratemap{17} = nan(size(pad_floorToP3_x)-1); coordmap_x{17} = pad_floorToP3_x; coordmap_y{17} = pad_floorToP3_y;
    ratemap{18} = nan(size(pad_floorToP4_x)-1); coordmap_x{18} = pad_floorToP4_x; coordmap_y{18} = pad_floorToP4_y;
    
    ratemap{19} = nan(size(pad_P1BToP1L_x)-1); coordmap_x{19} = pad_P1BToP1L_x; coordmap_y{19} = pad_P1BToP1L_y;
    ratemap{20} = nan(size(pad_P2BToP2L_x)-1); coordmap_x{20} = pad_P2BToP2L_x; coordmap_y{20} = pad_P2BToP2L_y;
    ratemap{21} = nan(size(pad_P3BToP3L_x)-1); coordmap_x{21} = pad_P3BToP3L_x; coordmap_y{21} = pad_P3BToP3L_y;
    ratemap{22} = nan(size(pad_P4BToP4L_x)-1); coordmap_x{22} = pad_P4BToP4L_x; coordmap_y{22} = pad_P4BToP4L_y;
    
    for ii = 1:7
        for row = 1:size(ratemap{ii},1)
            for col = 1:size(ratemap{ii},2)
                
                x = coordmap_x{ii}(row,col);
                y = coordmap_y{ii}(row,col);
                
                x = (x - 20) / 60;
                y = (y - 20) / 60;
                
                poiss_lambda = nan(svfield_num,1);
                for svf = 1:svfield_num
                    exponent = view_params((svf-1)*6 + 1)*x + view_params((svf-1)*6 + 2)*(x^2) + ...
                        view_params((svf-1)*6 + 3)*y + view_params((svf-1)*6 + 4)*(y^2) + ...
                        view_params((svf-1)*6 + 5)*x*y + view_params((svf-1)*6 + 6);
                    poiss_lambda(svf) = exp(exponent);
                end
                poiss_lambda = sum(poiss_lambda);
                
                ratemap{ii}(row,col) = poiss_lambda;
            end
        end
        ratemap{ii} = ratemap{ii} / tbin_size;
    end
    
    hold on
    
    for ii = 1:7
        coordmap_z = zeros(size(coordmap_x{ii}));
        coordmap_z(isnan(coordmap_x{ii})) = nan;
        
        surf(coordmap_x{ii}, coordmap_y{ii}, coordmap_z, ratemap{ii});
        alpha 1; shading flat;
        view(-35,20);
        colormap jet;
        colorbar;
    end

end

end
