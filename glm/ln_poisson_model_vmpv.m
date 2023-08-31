% this function is a slimmer version of the github reference version (since
% we only have view, place and head direction). linked below:
% https://github.com/GiocomoLab/ln-model-of-mec-neurons/blob/master/ln_poisson_model.m
%
% inserted 3 arguments to track place/head direction/view bins as well (for ease of use)

function [f, df, hessian] = ln_poisson_model_vmpv(param,data,modelType,np,nh,nv,beta_p,beta_h,beta_v)

    X = data{1}; % subset of A
    Y = data{2}; % number of spikes

    % compute the firing rate
    u = X * param;
    rate = exp(u);

    % roughness regularizer weight - note: these are tuned using the sum of f,
    % and thus have decreasing influence with increasing amounts of data
    b_pos = beta_p; b_hd = beta_h; b_view = beta_v;

    % start computing the Hessian
    rX = bsxfun(@times,rate,X);
    hessian_glm = rX'*X;

    %% find the parameters and compute their roughness penalties

    % initialize parameter-relevant variables
    J_pos = 0; J_pos_g = []; J_pos_h = [];
    J_hd = 0; J_hd_g = []; J_hd_h = [];
    J_view = 0; J_view_g = []; J_view_h = [];

    % find the parameters
    numPos = np; numHD = nh; numView = nv; % hardcoded: number of parameters
    [param_pos,param_hd,param_view] = find_param(param,modelType,numPos,numHD,numView);

    % compute the contribution for f, df, and the hessian
    if ~isempty(param_pos)
        [J_pos,J_pos_g,J_pos_h] = rough_penalty_pos(param_pos,b_pos);
    end
    
    if ~isempty(param_hd)
        [J_hd,J_hd_g,J_hd_h] = rough_penalty_hd(param_hd,b_hd);
    end

    if ~isempty(param_view)
        [J_view,J_view_g,J_view_h] = rough_penalty_spatialview(param_view,b_view);
    end

    %% compute f, the gradient, and the hessian

    f = sum(rate-Y.*u) + J_pos + J_hd + J_view;
    df = real(X' * (rate - Y) + [J_pos_g; J_hd_g; J_view_g]);
    hessian = hessian_glm + blkdiag(J_pos_h,J_hd_h,J_view_h);

end
    
    
%% smoothing functions called in the above script
function [J,J_g,J_h] = rough_penalty_pos(param,beta)

    numParam = numel(param);
    D1 = spdiags(ones(sqrt(numParam),1)*[-1 1],0:1,sqrt(numParam)-1,sqrt(numParam));
    DD1 = D1'*D1;
    M1 = kron(eye(sqrt(numParam)),DD1); M2 = kron(DD1,eye(sqrt(numParam)));
    M = (M1 + M2);

    J = beta*0.5*param'*M*param;
    J_g = beta*M*param;
    J_h = beta*M;

end

function [J,J_g,J_h] = rough_penalty_hd(param,beta)
    
    numParam = numel(param);
    D1 = spdiags(ones(numParam,1)*[-1 1],0:1,numParam-1,numParam);
    DD1 = D1'*D1;
    
    % to correct the smoothing across first and last bin
    DD1(1,:) = circshift(DD1(2,:),[0 -1]);
    DD1(end,:) = circshift(DD1(end-1,:),[0 1]);
    
    J = beta*0.5*param'*DD1*param;
    J_g = beta*DD1*param;
    J_h = beta*DD1;
    
end

function [J,J_g,J_h] = rough_penalty_spatialview(param,beta)

    % params increase across columns, then rows.
    % e.g.
    % 1 2 3 4 5
    % 6 7 8 9 10

    % floor and ceiling are 40 rows by 40 columns
    floor_width = 40;
    D1_floor = spdiags(ones(floor_width,1)*[-1 1],0:1,floor_width-1,floor_width);
    DD1_floor = D1_floor'*D1_floor;
    M1_floor = kron(eye(floor_width),DD1_floor); M2_floor = kron(DD1_floor,eye(floor_width));
    M_floor = (M1_floor + M2_floor);
    
    % walls are 8 rows by 40*4 columns
    % D1 for roughness across columns, D2 for roughness across rows
    wall_rows = 8; wall_cols = 40*4;
    D1_walls = spdiags(ones(wall_cols,1)*[-1 1],0:1,wall_cols-1,wall_cols);
    DD1_walls = D1_walls'*D1_walls;
    % to correct the smoothing across first and last wall edges
    DD1_walls(1,:) = circshift(DD1_walls(2,:),[0 -1]);
    DD1_walls(end,:) = circshift(DD1_walls(end-1,:),[0 1]);
    D2_walls = spdiags(ones(wall_rows,1)*[-1 1],0:1,wall_rows-1,wall_rows);
    DD2_walls = D2_walls'*D2_walls;
    M1_walls = kron(eye(wall_rows),DD1_walls); M2_walls = kron(DD2_walls,eye(wall_cols));
    M_walls = (M1_walls + M2_walls);
    
    % pillars are all 5 rows by 8*4 columns
    % D1 for roughness across columns, D2 for roughness across rows
    pillar_rows = 5; pillar_cols = 8*4;
    D1_pillars = spdiags(ones(pillar_cols,1)*[-1 1],0:1,pillar_cols-1,pillar_cols);
    DD1_pillars = D1_pillars'*D1_pillars;
    % to correct the smoothing across first and last wall edges
    DD1_pillars(1,:) = circshift(DD1_pillars(2,:),[0 -1]);
    DD1_pillars(end,:) = circshift(DD1_pillars(end-1,:),[0 1]);
    D2_pillars = spdiags(ones(pillar_rows,1)*[-1 1],0:1,pillar_rows-1,pillar_rows);
    DD2_pillars = D2_pillars'*D2_pillars;
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
        M(bin1, bin1) = M(bin1, bin1) + 1;
        M(bin2, bin2) = M(bin2, bin2) + 1;
        M(bin1, bin2) = M(bin1, bin2) - 1;
        M(bin2, bin1) = M(bin2, bin1) - 1;
        
        % ceiling and walls
        bin1 = ceiling_edge(i); bin2 = wall_top_edge(i);
        M(bin1, bin1) = M(bin1, bin1) + 1;
        M(bin2, bin2) = M(bin2, bin2) + 1;
        M(bin1, bin2) = M(bin1, bin2) - 1;
        M(bin2, bin1) = M(bin2, bin1) - 1;
    end
    
    % smoothing penalties between pillars and floor
    for i = 1:length(P1_BR_edge)
        % pillar 1
        bin1 = P1_BR_edge(i); bin2 = floor_P1_edge(i);
        M(bin1, bin1) = M(bin1, bin1) + 1;
        M(bin2, bin2) = M(bin2, bin2) + 1;
        M(bin1, bin2) = M(bin1, bin2) - 1;
        M(bin2, bin1) = M(bin2, bin1) - 1;
        
        % pillar 2
        bin1 = P2_BL_edge(i); bin2 = floor_P2_edge(i);
        M(bin1, bin1) = M(bin1, bin1) + 1;
        M(bin2, bin2) = M(bin2, bin2) + 1;
        M(bin1, bin2) = M(bin1, bin2) - 1;
        M(bin2, bin1) = M(bin2, bin1) - 1;
        
        % pillar 3
        bin1 = P3_TR_edge(i); bin2 = floor_P3_edge(i);
        M(bin1, bin1) = M(bin1, bin1) + 1;
        M(bin2, bin2) = M(bin2, bin2) + 1;
        M(bin1, bin2) = M(bin1, bin2) - 1;
        M(bin2, bin1) = M(bin2, bin1) - 1;
        
        % pillar 4
        bin1 = P4_TL_edge(i); bin2 = floor_P4_edge(i);
        M(bin1, bin1) = M(bin1, bin1) + 1;
        M(bin2, bin2) = M(bin2, bin2) + 1;
        M(bin1, bin2) = M(bin1, bin2) - 1;
        M(bin2, bin1) = M(bin2, bin1) - 1;
    end

    J = beta*0.5*param'*M*param;
    J_g = beta*M*param;
    J_h = beta*M;
    
end

%% function to find the right parameters given the model type
function [param_pos,param_hd,param_view] = find_param(param,modelType,numPos,numHD,numView)

    param_pos = []; param_hd = []; param_view = [];

    if all(modelType == [1 0 0])
        param_pos = param;
    elseif all(modelType == [0 1 0])
        param_hd = param;
    elseif all(modelType == [0 0 1])
        param_view = param;
    elseif all(modelType == [1 1 0])
        param_pos = param(1:numPos);
        param_hd = param(numPos+1:numPos+numHD);
    elseif all(modelType == [1 0 1])
        param_pos = param(1:numPos);
        param_view = param(numPos+1:numPos+numView);
    elseif all(modelType == [0 1 1])
        param_hd = param(1:numHD);
        param_view = param(numHD+1:numHD+numView);
    elseif all(modelType == [1 1 1])
        param_pos = param(1:numPos);
        param_hd = param(numPos+1:numPos+numHD);
        param_view = param(numPos+numHD+1:numPos+numHD+numView);
    end

end
