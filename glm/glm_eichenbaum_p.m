function glm_eichenbaum_p(glm_data, max_pfields, max_svfields, pillcon)
%
%   Reference: Macdonald, Lepage, Eden, and Eichenbaum, 2011. Hippocampal eeTime Cellsff Bridge the Gap in Memory for Discontiguous Events
%
%   Place field parameters are implemented the same as reference paper,
%   except for one extra parameter allowing for overall magnitude of field
%   to vary.
%   For fitting the spatialview fields, the view bins have been flattened
%   into a 2D plane. The pillars, walls, ceiling etc. are connected to each
%   other with padding to allow for fields spanning across these areas.
%
%   (SV PADDING HAS BEEN DISABLED FOR NOW)
%
%   Run this function in the cell directory, same folder as spiketrain.mat.
%
%   PARAMETERS:
%   glm_data - uses either glm_vmpvData or glm_genData.
%
%   max_pfields - maximum number of gaussian-shaped place fields to fit. Fitting
%   will be performed multiple times for 1:max_pfields and the best fit can
%   be decided by AIC.
%
%   max_svfields - maximum number of gaussian-shaped spatialview fields to fit.
%
%   pillcon - optional parameter, set to true to use constrained fitting
%   for the place or view model to disallow fitted means above a certain limit, in
%   the bins where the pillars should be. This was found to dampen
%   the overall intensity of the fitted fields, so it might not be a good
%   idea to use this.
%   

eb_results = struct;

bin_stc = glm_data.bin_stc;
tbin_size = glm_data.tbin_size;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mapping coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Code adapted from plotgridmap.m

padding = 5;

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

% Pillar areas, to set fitting constraints
P1a_x = repmat(24:31, 8, 1);
P1a_y = repmat([8:15]', 1, 8);

P2a_x = repmat(8:15, 8, 1);
P2a_y = repmat([8:15]', 1, 8);

P3a_x = repmat(24:31, 8, 1);
P3a_y = repmat([24:31]', 1, 8);

P4a_x = repmat(8:15, 8, 1);
P4a_y = repmat([24:31]', 1, 8);

pillar_a_x = [P1a_x(:); P2a_x(:); P3a_x(:); P4a_x(:)];
pillar_a_y = [P1a_y(:); P2a_y(:); P3a_y(:); P4a_y(:)];

% Constraints for spacing between meshes and at corners of meshes, for
% spatial view
% Used to stop fields from spanning from one mesh to the other, or
% "snaking" around corners.

pill_spacing = 20-5-padding; % typical y distance between pillar meshes

space_P1_P2_y = (sv_P1_y(1,1)+pill_spacing/2)*ones(1, 8*4 + 1);
space_P1_P2_x = [-120:-88];

space_P2_P3_y = (sv_P2_y(1,1)+pill_spacing/2)*ones(1, 8*4 + 1);
space_P2_P3_x = space_P1_P2_x;

space_P3_P4_y = (sv_P3_y(1,1)+pill_spacing/2)*ones(1, 8*4 + 1);
space_P3_P4_x = space_P1_P2_x;

space_main_pill_y = [0:88]'; % long vertical line separating main mesh from pillar meshes
space_main_pill_x = 70*ones(size(space_main_pill_y));

space_diag_BR_x = [40:40+39];
space_diag_BR_y = [sv_walls_y(end,end)-padding:-1:sv_walls_y(end,end)-padding-39];

space_diag_BL_x = [-1:-1:-1-39];
space_diag_BL_y = space_diag_BR_y;

space_diag_TR_x = space_diag_BR_x;
space_diag_TR_y = [sv_walls_y(1,1)+padding+1:sv_walls_y(1,1)+padding+1+39];

space_diag_TL_x = space_diag_BL_x;
space_diag_TL_y = space_diag_TR_y;

spacing_x = [space_P1_P2_x(:); space_P2_P3_x(:); space_P3_P4_x(:); space_main_pill_x(:); space_diag_BR_x(:); space_diag_BL_x(:); space_diag_TR_x(:); space_diag_TL_x(:)];
spacing_y = [space_P1_P2_y(:); space_P2_P3_y(:); space_P3_P4_y(:); space_main_pill_y(:); space_diag_BR_y(:); space_diag_BL_y(:); space_diag_TR_y(:); space_diag_TL_y(:)];

%%%
sv_con_x = [pillar_a_x; spacing_x];
sv_con_y = [pillar_a_y; spacing_y];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sample_data_indep = nan(size(bin_stc,1), 6);
sample_data_indep(:,1) = bin_stc(:,1);

sample_data_dep = nan(size(bin_stc,1),1);
sample_data_dep(:,1) = bin_stc(:,4);

invalid_rows = [];

for k = 1:size(bin_stc,1)
    
    bin_place = bin_stc(k,2);
    [row,col] = find(place_map == bin_place);
    sample_data_indep(k,2) = floor_x(row,col);
    sample_data_indep(k,3) = floor_y(row,col);
    
    bin_view = bin_stc(k,3);
    sample_data_indep(k,6) = bin_view;
    if any(3:3+1600-1 == bin_view) % floor
        [row,col] = find(floor == bin_view);
        sample_data_indep(k,4) = sv_floor_x(row,col);
        sample_data_indep(k,5) = sv_floor_y(row,col);
%         if find(link_floorToWallL_x==sample_data_indep(k,4) & link_floorToWallL_y==sample_data_indep(k,5))
%             sample_data_indep = [sample_data_indep; padded_view_data(sample_data_indep(k,:), link_floorToWallL_x, link_floorToWallL_y, pad_floorToWallL_x, pad_floorToWallL_y)];
%             sample_data_dep = [sample_data_dep; sample_data_dep(k)];
%         elseif find(link_floorToWallR_x==sample_data_indep(k,4) & link_floorToWallR_y==sample_data_indep(k,5))
%             sample_data_indep = [sample_data_indep; padded_view_data(sample_data_indep(k,:), link_floorToWallR_x, link_floorToWallR_y, pad_floorToWallR_x, pad_floorToWallR_y)];
%             sample_data_dep = [sample_data_dep; sample_data_dep(k)];
%         elseif find(link_floorToWallB_x==sample_data_indep(k,4) & link_floorToWallB_y==sample_data_indep(k,5))
%             sample_data_indep = [sample_data_indep; padded_view_data(sample_data_indep(k,:), link_floorToWallB_x, link_floorToWallB_y, pad_floorToWallB_x, pad_floorToWallB_y)];
%             sample_data_dep = [sample_data_dep; sample_data_dep(k)];
%         elseif find(link_floorToP1_x==sample_data_indep(k,4) & link_floorToP1_y==sample_data_indep(k,5))
%             sample_data_indep = [sample_data_indep; padded_view_data(sample_data_indep(k,:), link_floorToP1_x, link_floorToP1_y, pad_floorToP1_x, pad_floorToP1_y)];
%             sample_data_dep = [sample_data_dep; sample_data_dep(k)];
%         elseif find(link_floorToP2_x==sample_data_indep(k,4) & link_floorToP2_y==sample_data_indep(k,5))
%             sample_data_indep = [sample_data_indep; padded_view_data(sample_data_indep(k,:), link_floorToP2_x, link_floorToP2_y, pad_floorToP2_x, pad_floorToP2_y)];
%             sample_data_dep = [sample_data_dep; sample_data_dep(k)];
%         elseif find(link_floorToP3_x==sample_data_indep(k,4) & link_floorToP3_y==sample_data_indep(k,5))
%             sample_data_indep = [sample_data_indep; padded_view_data(sample_data_indep(k,:), link_floorToP3_x, link_floorToP3_y, pad_floorToP3_x, pad_floorToP3_y)];
%             sample_data_dep = [sample_data_dep; sample_data_dep(k)];
%         elseif find(link_floorToP4_x==sample_data_indep(k,4) & link_floorToP4_y==sample_data_indep(k,5))
%             sample_data_indep = [sample_data_indep; padded_view_data(sample_data_indep(k,:), link_floorToP4_x, link_floorToP4_y, pad_floorToP4_x, pad_floorToP4_y)];
%             sample_data_dep = [sample_data_dep; sample_data_dep(k)];
%         end
    elseif any(1603:1603+1600-1 == bin_view) % ceiling
        [row,col] = find(ceiling == bin_view);
        sample_data_indep(k,4) = sv_ceiling_x(row,col);
        sample_data_indep(k,5) = sv_ceiling_y(row,col);
%         if find(link_ceilingToWallL_x==sample_data_indep(k,4) & link_ceilingToWallL_y==sample_data_indep(k,5))
%             sample_data_indep = [sample_data_indep; padded_view_data(sample_data_indep(k,:), link_ceilingToWallL_x, link_ceilingToWallL_y, pad_ceilingToWallL_x, pad_ceilingToWallL_y)];
%             sample_data_dep = [sample_data_dep; sample_data_dep(k)];
%         elseif find(link_ceilingToWallR_x==sample_data_indep(k,4) & link_ceilingToWallR_y==sample_data_indep(k,5))
%             sample_data_indep = [sample_data_indep; padded_view_data(sample_data_indep(k,:), link_ceilingToWallR_x, link_ceilingToWallR_y, pad_ceilingToWallR_x, pad_ceilingToWallR_y)];
%             sample_data_dep = [sample_data_dep; sample_data_dep(k)];
%         elseif find(link_ceilingToWallB_x==sample_data_indep(k,4) & link_ceilingToWallB_y==sample_data_indep(k,5))
%             sample_data_indep = [sample_data_indep; padded_view_data(sample_data_indep(k,:), link_ceilingToWallB_x, link_ceilingToWallB_y, pad_ceilingToWallB_x, pad_ceilingToWallB_y)];
%             sample_data_dep = [sample_data_dep; sample_data_dep(k)];
%         end
    elseif any(3203:3203+1280-1 == bin_view) % walls
        [row,col] = find(walls == bin_view);
        sample_data_indep(k,4) = sv_walls_x(row,col);
        sample_data_indep(k,5) = sv_walls_y(row,col);
%         if find(link_wallBToWallL_x==sample_data_indep(k,4) & link_wallBToWallL_y==sample_data_indep(k,5))
%             sample_data_indep = [sample_data_indep; padded_view_data(sample_data_indep(k,:), link_wallBToWallL_x, link_wallBToWallL_y, pad_wallBToWallL_x, pad_wallBToWallL_y)];
%             sample_data_dep = [sample_data_dep; sample_data_dep(k)];
%         end
    elseif any(4483:4483+160-1 == bin_view) % P1_BR
        [row,col] = find(P1_BR == bin_view);
        sample_data_indep(k,4) = sv_P1_x(row,col);
        sample_data_indep(k,5) = sv_P1_y(row,col);
%         if find(link_P1BToP1L_x==sample_data_indep(k,4) & link_P1BToP1L_y==sample_data_indep(k,5))
%             sample_data_indep = [sample_data_indep; padded_view_data(sample_data_indep(k,:), link_P1BToP1L_x, link_P1BToP1L_y, pad_P1BToP1L_x, pad_P1BToP1L_y)];
%             sample_data_dep = [sample_data_dep; sample_data_dep(k)];
%         end
    elseif any(4643:4643+160-1 == bin_view) % P2_BL
        [row,col] = find(P2_BL == bin_view);
        sample_data_indep(k,4) = sv_P2_x(row,col);
        sample_data_indep(k,5) = sv_P2_y(row,col);
%         if find(link_P2BToP2L_x==sample_data_indep(k,4) & link_P2BToP2L_y==sample_data_indep(k,5))
%             sample_data_indep = [sample_data_indep; padded_view_data(sample_data_indep(k,:), link_P2BToP2L_x, link_P2BToP2L_y, pad_P2BToP2L_x, pad_P2BToP2L_y)];
%             sample_data_dep = [sample_data_dep; sample_data_dep(k)];
%         end
    elseif any(4803:4803+160-1 == bin_view)
        [row,col] = find(P3_TR == bin_view); % P3_TR
        sample_data_indep(k,4) = sv_P3_x(row,col);
        sample_data_indep(k,5) = sv_P3_y(row,col);
%         if find(link_P3BToP3L_x==sample_data_indep(k,4) & link_P3BToP3L_y==sample_data_indep(k,5))
%             sample_data_indep = [sample_data_indep; padded_view_data(sample_data_indep(k,:), link_P3BToP3L_x, link_P3BToP3L_y, pad_P3BToP3L_x, pad_P3BToP3L_y)];
%             sample_data_dep = [sample_data_dep; sample_data_dep(k)];
%         end
    elseif any(4963:4963+160-1 == bin_view) % P4_TL
        [row,col] = find(P4_TL == bin_view);
        sample_data_indep(k,4) = sv_P4_x(row,col);
        sample_data_indep(k,5) = sv_P4_y(row,col);
%         if find(link_P4BToP4L_x==sample_data_indep(k,4) & link_P4BToP4L_y==sample_data_indep(k,5))
%             sample_data_indep = [sample_data_indep; padded_view_data(sample_data_indep(k,:), link_P4BToP4L_x, link_P4BToP4L_y, pad_P4BToP4L_x, pad_P4BToP4L_y)];
%             sample_data_dep = [sample_data_dep; sample_data_dep(k)];
%         end
    else % cue or hint
        invalid_rows = [invalid_rows k];
    end
    
end

sample_data_indep(invalid_rows,:) = [];
sample_data_dep(invalid_rows,:) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEMONSTRATION PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A flattened firing rate map of the glm_vmpvData.
% To be used to compare with the raw sv rate maps. plot_rolloutsv can be
% used for comparison.

figure
spikes_dur = nan(5122,2);
for k = 1:size(sample_data_indep,1)
    bin_view = sample_data_indep(k,6);
    if isnan(spikes_dur(bin_view,2))
        spikes_dur(bin_view,1) = 0;
        spikes_dur(bin_view,2) = 0;
    end
    spikes_dur(bin_view,1) = spikes_dur(bin_view,1) + sample_data_dep(k,1);
    spikes_dur(bin_view,2) = spikes_dur(bin_view,2) + 1;
end
ratemap = spikes_dur(:,1) ./ spikes_dur(:,2);
ratemap = ratemap ./ tbin_size;


% Plot floor
surf(sv_floor_x, sv_floor_y, zeros(size(sv_floor_x)), flipud(reshape(ratemap(3:1600+3-1), 40, 40)'));
alpha 1; shading flat;
hold on;

% Plot ceiling and walls
surf(sv_ceiling_x, sv_ceiling_y, zeros(size(sv_ceiling_y)), flipud(reshape(ratemap(1603:1603+1600-1), 40, 40)'));
alpha 1; shading flat;
surf(sv_walls_x, sv_walls_y, zeros(size(sv_walls_x)), flipud(reshape(ratemap(3203:3203+1280-1), 40*4, 8)'));      
alpha 1; shading flat;

% Plot pillars
surf(sv_P1_x, sv_P1_y, zeros(size(sv_P1_x)), flipud(reshape(ratemap(4483:4483+160-1), 8*4, 5)'));
alpha 1; shading flat;
surf(sv_P2_x, sv_P2_y, zeros(size(sv_P1_x)), flipud(reshape(ratemap(4643:4643+160-1), 8*4, 5)'));
alpha 1; shading flat;
surf(sv_P3_x, sv_P3_y, zeros(size(sv_P1_x)), flipud(reshape(ratemap(4803:4803+160-1), 8*4, 5)'));
alpha 1; shading flat;
surf(sv_P4_x, sv_P4_y, zeros(size(sv_P1_x)), flipud(reshape(ratemap(4963:4963+160-1), 8*4, 5)'));
alpha 1; shading flat; 
view(-35,20);
colormap jet;
colorbar;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%options_unc = optimoptions('fminunc','MaxFunctionEvaluations',10000,'FiniteDifferenceType','central','FunctionTolerance',1e-9,'StepTolerance',1e-10);
options_unc = optimoptions('fminunc','SpecifyObjectiveGradient',true,'FiniteDifferenceType','central','FunctionTolerance',1e-9,'StepTolerance',1e-10);
options_con = optimoptions('fmincon','MaxFunctionEvaluations',10000,'SpecifyObjectiveGradient',true,'FiniteDifferenceType','central','FunctionTolerance',1e-9,'StepTolerance',1e-12);
%options_con = optimoptions('fmincon','MaxFunctionEvaluations',10000,'FiniteDifferenceType','central','FunctionTolerance',1e-9,'StepTolerance',1e-10);

con_ub = 0.30 * tbin_size; % max allowed rate inside pillars, for constrained case

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Place model fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

place_params = cell(max_pfields, 1);
place_llh = cell(max_pfields, 1);

for pf = 1:max_pfields
    init_params = rand(pf*6,1) - 0.5;
    if nargin > 3 && pillcon
        A = []; B = []; Aeq = []; Beq = []; lb = []; ub = [];
        [temp_params, temp_llh] = fmincon(@(x)compute_llh_neg(x, pf, 0, sample_data_indep, sample_data_dep, 'place'), init_params, A, B, Aeq, Beq, lb, ub, @(x)pillar_cineq(x, pf, 0, pillar_a_x, pillar_a_y, con_ub, 'place'), options_con);
    else
        [temp_params, temp_llh] = fminunc(@(x)compute_llh_neg(x, pf, 0, sample_data_indep, sample_data_dep, 'place'), init_params, options_unc);
    end
    place_params{pf,1} = temp_params;
    place_llh{pf,1} = -temp_llh;
end

eb_results.params_place = place_params;
eb_results.llh_place = place_llh;
disp('place model done')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% View model fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

view_params = cell(1, max_svfields);
view_llh = cell(1, max_svfields);

for svf = 1:max_svfields
    init_params = rand(svf*6,1) - 0.5;
    if nargin > 3 && pillcon
        A = []; B = []; Aeq = []; Beq = []; lb = []; ub = [];
        [temp_params, temp_llh] = fmincon(@(x)compute_llh_neg(x, 0, svf, sample_data_indep, sample_data_dep, 'view'), init_params, A, B, Aeq, Beq, lb, ub, @(x)pillar_cineq(x, 0, svf, sv_con_x, sv_con_y, con_ub, 'view'), options_con);
    else
        [temp_params, temp_llh] = fminunc(@(x)compute_llh_neg(x, 0, svf, sample_data_indep, sample_data_dep, 'view'), init_params, options_unc);
    end
    view_params{1, svf} = temp_params;
    view_llh{1, svf} = -temp_llh;
end

eb_results.params_view = view_params;
eb_results.llh_view = view_llh;
disp('view model done')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Joint model fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

joint_params = cell(max_pfields, max_svfields);
joint_llh = cell(max_pfields, max_svfields);

for pf = 1:max_pfields
    for svf = 1:max_svfields
        init_params = rand(pf*6+svf*6,1) - 0.5;
        if nargin > 3 && pillcon
            A = []; B = []; Aeq = []; Beq = []; lb = []; ub = [];
            %[temp_params, temp_llh] = fmincon(@(x)compute_llh_neg(x, pf, svf, sample_data_indep, sample_data_dep, 'both'), init_params, A, B, Aeq, Beq, lb, ub, @(x)pillar_cineq(x, pf, svf, sv_con_x, sv_con_y, con_ub, 'both'), options_con);
            [temp_params, temp_llh] = fminunc(@(x)compute_llh_neg(x, pf, svf, sample_data_indep, sample_data_dep, 'both'), init_params, options_unc);
        else
            [temp_params, temp_llh] = fminunc(@(x)compute_llh_neg(x, pf, svf, sample_data_indep, sample_data_dep, 'both'), init_params, options_unc);
        end
        joint_params{pf, svf} = temp_params;
        joint_llh{pf, svf} = -temp_llh;
    end
end

eb_results.params_joint = joint_params;
eb_results.llh_joint = joint_llh;
disp('joint model done')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W_place = cell(max_pfields, max_svfields);
W_view = cell(max_pfields, max_svfields);
STIC_placeoverview = cell(max_pfields, max_svfields);
AIC_place = cell(max_pfields, 1);
AIC_view = cell(1, max_svfields);
AIC_joint = cell(max_pfields, max_svfields);


for pf = 1:max_pfields
    for svf = 1:max_svfields
        W_place{pf, svf} = -2*(place_llh{pf,1} - joint_llh{pf, svf});
        W_view{pf, svf} = -2*(view_llh{1,svf} - joint_llh{pf, svf});
        STIC_placeoverview{pf,svf} = place_llh{pf,1} - view_llh{1,svf};
        AIC_joint{pf, svf} = 2*(6*pf + 10*svf) - 2*joint_llh{pf, svf};
    end
end

for pf = 1:max_pfields
    AIC_place{pf, 1} = 2*(6*pf) - 2*place_llh{pf,1};
end
for svf = 1:max_svfields
    AIC_view{1, svf} = 2*(10*svf) - 2*view_llh{1,svf};
end

eb_results.W_place = W_place;
eb_results.W_view = W_view;
eb_results.STIC_placeoverview = STIC_placeoverview;
eb_results.tbin_size = glm_data.tbin_size;
eb_results.padding = padding;
eb_results.ThresVel = glm_data.ThresVel;
eb_results.UseMinObs = glm_data.UseMinObs;

eb_results.AIC_place = AIC_place;
eb_results.AIC_view = AIC_view;
eb_results.AIC_joint = AIC_joint;

if nargin > 3
    eb_results.pillcon = pillcon;
end


save('glm_eb_p_results_1_5_nopad_1ms_unc.mat','eb_results','-v7.3');

end

function [f, g] = compute_llh_neg(input_params, pfield_num, svfield_num, data_indep, data_dep, model)
%

llh_entries = nan(size(data_indep,1), 1);
gradient = nan(size(data_indep,1), length(input_params));

switch model
    case 'place' % 6 parameters
        for k = 1:size(llh_entries,1)
            place_x = (data_indep(k,2) - 20) / 60;
            place_y = (data_indep(k,3) - 20) / 60;
            
            poiss_lambda = nan(pfield_num,1);
            for pf = 1:pfield_num
                exponent = input_params((pf-1)*6 + 1)*place_x + input_params((pf-1)*6 + 2)*(place_x^2) + ...
                    input_params((pf-1)*6 + 3)*place_y + input_params((pf-1)*6 + 4)*(place_y^2) + ...
                    input_params((pf-1)*6 + 5)*place_x*place_y + input_params((pf-1)*6 + 6);
                if exp(exponent) == 0
                    poiss_lambda(pf) = 1e-14;
                else
                    poiss_lambda(pf) = exp(exponent);
                end
                if nargout > 1 % gradient demanded
                    gradient(k, (pf-1)*6 + 1) = poiss_lambda(pf) * place_x;
                    gradient(k, (pf-1)*6 + 2) = poiss_lambda(pf) * place_x^2;
                    gradient(k, (pf-1)*6 + 3) = poiss_lambda(pf) * place_y;
                    gradient(k, (pf-1)*6 + 4) = poiss_lambda(pf) * place_y^2;
                    gradient(k, (pf-1)*6 + 5) = poiss_lambda(pf) * place_x*place_y;
                    gradient(k, (pf-1)*6 + 6) = poiss_lambda(pf);
                end
            end
            poiss_lambda = sum(poiss_lambda);
            if nargout > 1
                gradient(k, :) = ((data_dep(k)/poiss_lambda) - 1) * gradient(k, :);
            end
            llh_entries(k) = data_dep(k)*log(poiss_lambda) - poiss_lambda - log(factorial(data_dep(k)));
        end
    case 'view' % 6 parameters
        for k = 1:size(llh_entries,1)
            view_x = (data_indep(k,4) - 20) / 60;
            view_y = (data_indep(k,5) - 20) / 60;
            
            poiss_lambda = nan(svfield_num,1);
            for svf = 1:svfield_num
                exponent = input_params((svf-1)*6 + 1)*view_x + input_params((svf-1)*6 + 2)*(view_x^2) + ...
                    input_params((svf-1)*6 + 3)*view_y + input_params((svf-1)*6 + 4)*(view_y^2) + ...
                    input_params((svf-1)*6 + 5)*view_x*view_y + input_params((svf-1)*6 + 6);
                if exp(exponent) == 0
                    poiss_lambda(svf) = 1e-14;
                else
                    poiss_lambda(svf) = exp(exponent);
                end
                if nargout > 1 % gradient demanded
                    gradient(k, (svf-1)*6 + 1) = poiss_lambda(svf) * view_x;
                    gradient(k, (svf-1)*6 + 2) = poiss_lambda(svf) * view_x^2;
                    gradient(k, (svf-1)*6 + 3) = poiss_lambda(svf) * view_y;
                    gradient(k, (svf-1)*6 + 4) = poiss_lambda(svf) * view_y^2;
                    gradient(k, (svf-1)*6 + 5) = poiss_lambda(svf) * view_x*view_y;
                    gradient(k, (svf-1)*6 + 6) = poiss_lambda(svf);
                end
            end
            poiss_lambda = sum(poiss_lambda);
            if nargout > 1
                gradient(k, :) = ((data_dep(k)/poiss_lambda) - 1) * gradient(k, :);
            end
            llh_entries(k) = data_dep(k)*log(poiss_lambda) - poiss_lambda - log(factorial(data_dep(k)));
        end
    case 'both' % 6+6 parameters
        for k = 1:size(llh_entries,1)
            place_x = (data_indep(k,2) - 20) / 60;
            place_y = (data_indep(k,3) - 20) / 60;
            view_x = (data_indep(k,4) - 20) / 60;
            view_y = (data_indep(k,5) - 20) / 60;
            
            place_params = input_params(1:6*pfield_num);
            view_params = input_params(6*pfield_num+1:end);
            
            poiss_lambda = nan(pfield_num+svfield_num,1);
            
            for pf = 1:pfield_num
                exponent = place_params((pf-1)*6 + 1)*place_x + place_params((pf-1)*6 + 2)*(place_x^2) + ...
                    place_params((pf-1)*6 + 3)*place_y + place_params((pf-1)*6 + 4)*(place_y^2) + ...
                    place_params((pf-1)*6 + 5)*place_x*place_y + place_params((pf-1)*6 + 6);
                if exp(exponent) == 0
                    poiss_lambda(pf) = 1e-14;
                else
                    poiss_lambda(pf) = exp(exponent);
                end
                if nargout > 1 % gradient demanded
                    gradient(k, (pf-1)*6 + 1) = poiss_lambda(pf) * place_x;
                    gradient(k, (pf-1)*6 + 2) = poiss_lambda(pf) * place_x^2;
                    gradient(k, (pf-1)*6 + 3) = poiss_lambda(pf) * place_y;
                    gradient(k, (pf-1)*6 + 4) = poiss_lambda(pf) * place_y^2;
                    gradient(k, (pf-1)*6 + 5) = poiss_lambda(pf) * place_x*place_y;
                    gradient(k, (pf-1)*6 + 6) = poiss_lambda(pf);
                end
            end
            
            for svf = 1:svfield_num
                exponent = exponent + ...
                    view_params((svf-1)*6 + 1)*view_x + view_params((svf-1)*6 + 2)*(view_x^2) + ...
                    view_params((svf-1)*6 + 3)*view_y + view_params((svf-1)*6 + 4)*(view_y^2) + ...
                    view_params((svf-1)*6 + 5)*view_x*view_y + view_params((svf-1)*6 + 6);
                if exp(exponent) == 0
                    poiss_lambda(pfield_num + svf) = 1e-14;
                else
                    poiss_lambda(pfield_num + svf) = exp(exponent);
                end
                if nargout > 1
                    gradient(k, 6*pfield_num + (svf-1)*6 + 1) = poiss_lambda(pfield_num + svf) * view_x;
                    gradient(k, 6*pfield_num + (svf-1)*6 + 2) = poiss_lambda(pfield_num + svf) * view_x^2;
                    gradient(k, 6*pfield_num + (svf-1)*6 + 3) = poiss_lambda(pfield_num + svf) * view_y;
                    gradient(k, 6*pfield_num + (svf-1)*6 + 4) = poiss_lambda(pfield_num + svf) * view_y^2;
                    gradient(k, 6*pfield_num + (svf-1)*6 + 5) = poiss_lambda(pfield_num + svf) * view_x*view_y;
                    gradient(k, 6*pfield_num + (svf-1)*6 + 6) = poiss_lambda(pfield_num + svf);
                end
            end
            
            % intermediate step
            if nargout > 1
                gradient(k, 1:6*pfield_num) = sum(poiss_lambda(pfield_num+1:end)) * gradient(k, 1:6*pfield_num);
                gradient(k, 6*pfield_num+1:end) = sum(poiss_lambda(1:pfield_num)) * gradient(k, 6*pfield_num+1:end);
            end
            
            poiss_lambda = sum(poiss_lambda(1:pfield_num)) * sum(poiss_lambda(pfield_num+1:end)); % place * view
            if nargout > 1
                gradient(k, :) = ((data_dep(k)/poiss_lambda) - 1) * gradient(k, :);
            end
            llh_entries(k) = data_dep(k)*log(poiss_lambda) - poiss_lambda - log(factorial(data_dep(k)));
        end
end

f = -sum(llh_entries);
if nargout > 1
    g = -sum(gradient, 1);
end

end

function [c, ceq] = pillar_cineq(input_params, pfield_num, svfield_num, x, y, upper_bound, model)

scaled_x = (x - 20) / 60;
scaled_y = (y - 20) / 60;

switch model
    case 'place' % 6 parameters
        poiss_lambda = nan(pfield_num,length(x));
        for pf = 1:pfield_num
            exponent = input_params((pf-1)*6 + 1)*scaled_x + input_params((pf-1)*6 + 2)*(scaled_x.^2) + ...
                input_params((pf-1)*6 + 3)*scaled_y + input_params((pf-1)*6 + 4)*(scaled_y.^2) + ...
                input_params((pf-1)*6 + 5)*scaled_x.*scaled_y + input_params((pf-1)*6 + 6);
            if exp(exponent) == 0
                poiss_lambda(pf,:) = 1e-14;
            else
                poiss_lambda(pf,:) = exp(exponent);
            end
        end
        poiss_lambda = sum(poiss_lambda,1);
    case 'view' % 6 parameters
        poiss_lambda = nan(svfield_num,length(x));
        for svf = 1:svfield_num
            exponent = input_params((svf-1)*6 + 1)*scaled_x + input_params((svf-1)*6 + 2)*(scaled_x.^2) + ...
                input_params((svf-1)*6 + 3)*scaled_y + input_params((svf-1)*6 + 4)*(scaled_y.^2) + ...
                input_params((svf-1)*6 + 5)*scaled_x.*scaled_y + input_params((svf-1)*6 + 6);
            if exp(exponent) == 0
                poiss_lambda(svf,:) = 1e-14;
            else
                poiss_lambda(svf,:) = exp(exponent);
            end
        end
        poiss_lambda = sum(poiss_lambda,1);
    case 'both' % 6+6 parameters
        error('Constraints not yet implemented for joint model.')
end


c = poiss_lambda - upper_bound;
ceq = [];

end

function output_row = padded_view_data(input_row, link_x, link_y, pad_x, pad_y)

[row, col] = find(link_x==input_row(4) & link_y==input_row(5));
output_row = input_row;
output_row(4) = pad_x(row,col);
output_row(5) = pad_y(row,col);

end