% this function is a slimmer version of the github reference version (since
% we only have view and place). linked below:
% https://github.com/GiocomoLab/ln-model-of-mec-neurons/blob/master/ln_poisson_model.m
%
% inserted 2 arguments to track place/view bins as well (for ease of use)

function [f, df, hessian] = ln_poisson_model_vmpv(param,data,modelType,np,nv,beta_p,beta_v)

    X = data{1}; % subset of A
    Y = data{2}; % number of spikes

    % compute the firing rate
    u = X * param;
    rate = exp(u);

    % roughness regularizer weight - note: these are tuned using the sum of f,
    % and thus have decreasing influence with increasing amounts of data
    b_pos = beta_p; b_view = beta_v;

    % start computing the Hessian
    rX = bsxfun(@times,rate,X);
    hessian_glm = rX'*X;

    %% find the parameters and compute their roughness penalties

    % initialize parameter-relevant variables
    J_pos = 0; J_pos_g = []; J_pos_h = [];
    J_view = 0; J_view_g = []; J_view_h = [];

    % find the parameters
    numPos = np; numView = nv; % hardcoded: number of parameters
    [param_pos,param_view] = find_param(param,modelType,numPos,numView);

    % compute the contribution for f, df, and the hessian
    if ~isempty(param_pos)
        [J_pos,J_pos_g,J_pos_h] = rough_penalty_pos(param_pos,b_pos);
    end

    if ~isempty(param_view)
        [J_view,J_view_g,J_view_h] = rough_penalty_spatialview(param_view,b_view);
    end

    %% compute f, the gradient, and the hessian

    f = sum(rate-Y.*u) + J_pos + J_view;
    df = real(X' * (rate - Y) + [J_pos_g; J_view_g]);
    hessian = hessian_glm + blkdiag(J_pos_h,J_view_h);

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

function [J,J_g,J_h] = rough_penalty_spatialview(param,beta)

    % params increase across columns, then rows.
    % e.g.
    % 1 2 3 4 5
    % 6 7 8 9 10

    param_floor = param(3:3+1600-1);
    numParam_floor = numel(param_floor);
    D1_floor = spdiags(ones(sqrt(numParam_floor),1)*[-1 1],0:1,sqrt(numParam_floor)-1,sqrt(numParam_floor));
    DD1_floor = D1_floor'*D1_floor;
    M1_floor = kron(eye(sqrt(numParam_floor)),DD1_floor); M2_floor = kron(DD1_floor,eye(sqrt(numParam_floor)));
    M_floor = (M1_floor + M2_floor);
    J_floor = beta*0.5*param_floor'*M_floor*param_floor;
    J_g_floor = beta*M_floor*param_floor;
    J_h_floor = beta*M_floor;
    
    param_ceiling = param(1603:1603+1600-1);
    numParam_ceiling = numel(param_ceiling);
    D1_ceiling = spdiags(ones(sqrt(numParam_ceiling),1)*[-1 1],0:1,sqrt(numParam_ceiling)-1,sqrt(numParam_ceiling));
    DD1_ceiling = D1_ceiling'*D1_ceiling;
    M1_ceiling = kron(eye(sqrt(numParam_ceiling)),DD1_ceiling); M2_ceiling = kron(DD1_ceiling,eye(sqrt(numParam_ceiling)));
    M_ceiling = (M1_ceiling + M2_ceiling);
    J_ceiling = beta*0.5*param_ceiling'*M_ceiling*param_ceiling;
    J_g_ceiling = beta*M_ceiling*param_ceiling;
    J_h_ceiling = beta*M_ceiling;   
    
    % walls are 8 rows by 40*4 columns
    param_walls = param(3203:3203+1280-1);
    % D1 for roughness across columns, D2 for roughness across rows
    D1_walls = spdiags(ones(40*4,1)*[-1 1],0:1,(40*4)-1,40*4);
    DD1_walls = D1_walls'*D1_walls;
    D2_walls = spdiags(ones(8,1)*[-1 1],0:1,8-1,8);
    DD2_walls = D2_walls'*D2_walls;
    M1_walls = kron(eye(8),DD1_walls); M2_walls = kron(DD2_walls,eye(40*4));
    M_walls = (M1_walls + M2_walls);
    J_walls = beta*0.5*param_walls'*M_walls*param_walls;
    J_g_walls = beta*M_walls*param_walls;
    J_h_walls = beta*M_walls;
    
    % pillars are all 5 rows by 8*4 columns
    
    param_P1_BR = param(4483:4483+160-1);
    D1_P1_BR = spdiags(ones(8*4,1)*[-1 1],0:1,(8*4)-1,8*4);
    DD1_P1_BR = D1_P1_BR'*D1_P1_BR;
    D2_P1_BR = spdiags(ones(5,1)*[-1 1],0:1,5-1,5);
    DD2_P1_BR = D2_P1_BR'*D2_P1_BR;
    M1_P1_BR = kron(eye(5),DD1_P1_BR); M2_P1_BR = kron(DD2_P1_BR,eye(8*4));
    M_P1_BR = (M1_P1_BR + M2_P1_BR);
    J_P1_BR = beta*0.5*param_P1_BR'*M_P1_BR*param_P1_BR;
    J_g_P1_BR = beta*M_P1_BR*param_P1_BR;
    J_h_P1_BR = beta*M_P1_BR;
    
    param_P2_BL = param(4643:4643+160-1);
    D1_P2_BL = spdiags(ones(8*4,1)*[-1 1],0:1,(8*4)-1,8*4);
    DD1_P2_BL = D1_P2_BL'*D1_P2_BL;
    D2_P2_BL = spdiags(ones(5,1)*[-1 1],0:1,5-1,5);
    DD2_P2_BL = D2_P2_BL'*D2_P2_BL;
    M1_P2_BL = kron(eye(5),DD1_P2_BL); M2_P2_BL = kron(DD2_P2_BL,eye(8*4));
    M_P2_BL = (M1_P2_BL + M2_P2_BL);
    J_P2_BL = beta*0.5*param_P2_BL'*M_P1_BR*param_P2_BL;
    J_g_P2_BL = beta*M_P2_BL*param_P2_BL;
    J_h_P2_BL = beta*M_P2_BL;
    
    param_P3_TR = param(4803:4803+160-1);
    D1_P3_TR = spdiags(ones(8*4,1)*[-1 1],0:1,(8*4)-1,8*4);
    DD1_P3_TR = D1_P3_TR'*D1_P3_TR;
    D2_P3_TR = spdiags(ones(5,1)*[-1 1],0:1,5-1,5);
    DD2_P3_TR = D2_P3_TR'*D2_P3_TR;
    M1_P3_TR = kron(eye(5),DD1_P3_TR); M2_P3_TR = kron(DD2_P3_TR,eye(8*4));
    M_P3_TR = (M1_P3_TR + M2_P3_TR);
    J_P3_TR = beta*0.5*param_P3_TR'*M_P1_BR*param_P3_TR;
    J_g_P3_TR = beta*M_P3_TR*param_P3_TR;
    J_h_P3_TR = beta*M_P3_TR;
    
    param_P4_TL = param(4963:4963+160-1);
    D1_P4_TL = spdiags(ones(8*4,1)*[-1 1],0:1,(8*4)-1,8*4);
    DD1_P4_TL = D1_P4_TL'*D1_P4_TL;
    D2_P4_TL = spdiags(ones(5,1)*[-1 1],0:1,5-1,5);
    DD2_P4_TL = D2_P4_TL'*D2_P4_TL;
    M1_P4_TL = kron(eye(5),DD1_P4_TL); M2_P4_TL = kron(DD2_P4_TL,eye(8*4));
    M_P4_TL = (M1_P4_TL + M2_P4_TL);
    J_P4_TL = beta*0.5*param_P4_TL'*M_P1_BR*param_P4_TL;
    J_g_P4_TL = beta*M_P4_TL*param_P4_TL;
    J_h_P4_TL = beta*M_P4_TL;
    
    J = J_floor + J_ceiling + J_walls + J_P1_BR + J_P2_BL + J_P3_TR + J_P4_TL;
    J_g = [0; 0; J_g_floor; J_g_ceiling; J_g_walls; J_g_P1_BR; J_g_P2_BL; J_g_P3_TR; J_g_P4_TL];
    J_h = blkdiag(zeros(2), J_h_floor, J_h_ceiling, J_h_walls, J_h_P1_BR, J_h_P2_BL, J_h_P3_TR, J_h_P4_TL);
    
end

%% function to find the right parameters given the model type
function [param_pos,param_view] = find_param(param,modelType,numPos,numView)

    param_pos = []; param_view = [];

    if all(modelType == [1 0])
        param_pos = param;
    elseif all(modelType == [0 1])
        param_view = param;
    elseif all(modelType == [1 1])
        param_pos = param(1:numPos);
        param_view = param(numPos+1:numPos+numView);
    end

end
