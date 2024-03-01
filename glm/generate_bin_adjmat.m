function M = generate_bin_adjmat(saveFile)
    
    if ~exist('saveFile', 'var')
        saveFile = false;
    end

    %%% Adjacency matrix for all view bins in the maze %%%
    % floor and ceiling are 40 rows by 40 columns
    floor_width = 40;
    DD1_floor = spdiags(ones(floor_width,1)*[1 1],[-1,1],floor_width,floor_width);
    M1_floor = kron(eye(floor_width),DD1_floor); M2_floor = kron(DD1_floor,eye(floor_width));
    M_floor = (M1_floor + M2_floor);
    
    % walls are 8 rows by 40*4 columns
    % D1 for roughness across columns, D2 for roughness across rows
    wall_rows = 8; wall_cols = 40*4;
    DD1_walls = spdiags(ones(wall_cols,1)*[1 1],[-1,1],wall_cols,wall_cols);
    % to correct the smoothing across first and last wall edges
    DD1_walls(1,:) = circshift(DD1_walls(2,:),[0 -1]);
    DD1_walls(end,:) = circshift(DD1_walls(end-1,:),[0 1]);
    DD2_walls = spdiags(ones(wall_rows,1)*[1 1],[-1,1],wall_rows,wall_rows);
    M1_walls = kron(eye(wall_rows),DD1_walls); M2_walls = kron(DD2_walls,eye(wall_cols));
    M_walls = (M1_walls + M2_walls);
    
    % pillars are all 5 rows by 8*4 columns
    % D1 for roughness across columns, D2 for roughness across rows
    pillar_rows = 5; pillar_cols = 8*4;
    DD1_pillars = spdiags(ones(pillar_cols,1)*[1 1],[-1,1],pillar_cols,pillar_cols);
    % to correct the smoothing across first and last wall edges
    DD1_pillars(1,:) = circshift(DD1_pillars(2,:),[0 -1]);
    DD1_pillars(end,:) = circshift(DD1_pillars(end-1,:),[0 1]);
    DD2_pillars = spdiags(ones(pillar_rows,1)*[1 1],[-1,1],pillar_rows,pillar_rows);
    M1_pillars = kron(eye(pillar_rows),DD1_pillars); M2_pillars = kron(DD2_pillars,eye(pillar_cols));
    M_pillars = (M1_pillars + M2_pillars);
    
    % combined smoothing penalty matrix
    M = blkdiag(zeros(2), M_floor, M_floor, M_walls, M_pillars, M_pillars, M_pillars, M_pillars);
    
    % get edges of floor, ceiling, walls and pillars
    floor_bins = reshape(1:floor_width^2, floor_width, floor_width) + 2;
    floor_edge = [floor_bins(1,:), floor_bins(:,end)', floor_bins(end,end:-1:1), floor_bins(end:-1:1,1)'];
    ceiling_edge = floor_edge + floor_width^2;
    wall_bottom_edge = (1:wall_cols) + 2*floor_width^2 + 2;
    wall_top_edge = wall_bottom_edge + (wall_rows-1)*wall_cols;
    P1_BR_edge = (1:pillar_cols) + wall_rows*wall_cols + 2*floor_width^2 + 2;
    P2_BL_edge = P1_BR_edge + pillar_rows*pillar_cols;
    P3_TR_edge = P2_BL_edge + pillar_rows*pillar_cols;
    P4_TL_edge = P3_TR_edge + pillar_rows*pillar_cols;
    
    % get edges of floor around pillars
    pillar_width = pillar_cols/4;
    P1_corner = [24, 8]; P2_corner = [8, 8]; P3_corner = [24, 24]; P4_corner = [8, 24];
    floor_P1_edge = floor_bins(P1_corner(1):P1_corner(1)+pillar_width+1, P1_corner(2):P1_corner(2)+pillar_width+1);
    floor_P2_edge = floor_bins(P2_corner(1):P2_corner(1)+pillar_width+1, P2_corner(2):P2_corner(2)+pillar_width+1);
    floor_P3_edge = floor_bins(P3_corner(1):P3_corner(1)+pillar_width+1, P3_corner(2):P3_corner(2)+pillar_width+1);
    floor_P4_edge = floor_bins(P4_corner(1):P4_corner(1)+pillar_width+1, P4_corner(2):P4_corner(2)+pillar_width+1);
    floor_P1_edge = [floor_P1_edge(1,2:end-1), floor_P1_edge(2:end-1,end)', floor_P1_edge(end,end-1:-1:2), floor_P1_edge(end-1:-1:2,1)'];
    floor_P2_edge = [floor_P2_edge(1,2:end-1), floor_P2_edge(2:end-1,end)', floor_P2_edge(end,end-1:-1:2), floor_P2_edge(end-1:-1:2,1)'];
    floor_P3_edge = [floor_P3_edge(1,2:end-1), floor_P3_edge(2:end-1,end)', floor_P3_edge(end,end-1:-1:2), floor_P3_edge(end-1:-1:2,1)'];
    floor_P4_edge = [floor_P4_edge(1,2:end-1), floor_P4_edge(2:end-1,end)', floor_P4_edge(end,end-1:-1:2), floor_P4_edge(end-1:-1:2,1)'];
    
    % smoothing penalties between floor/ceiling and walls
    for i = 1:length(floor_edge)
        % floor and walls
        bin1 = floor_edge(i); bin2 = wall_bottom_edge(i);
        M(bin1, bin2) = 1;
        M(bin2, bin1) = 1;
        
        % ceiling and walls
        bin1 = ceiling_edge(i); bin2 = wall_top_edge(i);
        M(bin1, bin2) = 1;
        M(bin2, bin1) = 1;
    end
    
    % smoothing penalties between pillars and floor
    for i = 1:length(P1_BR_edge)
        % pillar 1
        bin1 = P1_BR_edge(i); bin2 = floor_P1_edge(i);
        M(bin1, bin2) = 1;
        M(bin2, bin1) = 1;
        
        % pillar 2
        bin1 = P2_BL_edge(i); bin2 = floor_P2_edge(i);
        M(bin1, bin2) = 1;
        M(bin2, bin1) = 1;
        
        % pillar 3
        bin1 = P3_TR_edge(i); bin2 = floor_P3_edge(i);
        M(bin1, bin2) = 1;
        M(bin2, bin1) = 1;
        
        % pillar 4
        bin1 = P4_TL_edge(i); bin2 = floor_P4_edge(i);
        M(bin1, bin2) = 1;
        M(bin2, bin1) = 1;
    end
    
    % mark out view bins on the floor under pillars, and then remove
    % smoothing penalty for those view bins
    under_pillars = nan(256, 1); j = 1;
    for i = [9:16, 25:32]
        under_pillars(j:j+15) = 40*(i-1) + [9:16, 25:32];
        j = j+16;
    end
    under_pillars = under_pillars + 2;
    for i = 1:length(under_pillars)
        bin = under_pillars(i);
        M(bin, :) = 0;
        M(:, bin) = 0;
    end
    
    if saveFile
        save('doubletee_binmap.mat', 'M');
    end
    
end