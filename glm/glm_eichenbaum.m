function glm_eichenbaum(glm_data, max_pfields, max_svfields, pillcon)
%
%   Reference: Macdonald, Lepage, Eden, and Eichenbaum, 2011. Hippocampal eeTime Cellsff Bridge the Gap in Memory for Discontiguous Events
%
%   Place field parameters are implemented the same as reference paper,
%   except for one extra parameter allowing for overall magnitude of field
%   to vary.
%   Spatialview field parameters for this function are a simple extension
%   of the place field version to a 3D gaussian field.
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
%   for the place model to disallow fitted means above a certain limit, in
%   the place bins where the pillars should be. This was found to dampen
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
    if any(3:3+1600-1 == bin_view) % floor
        [row,col] = find(floor == bin_view);
        sample_data_indep(k,4) = floor_x(row,col);
        sample_data_indep(k,5) = floor_y(row,col);
        sample_data_indep(k,6) = floor_z(row,col);
    elseif any(1603:1603+1600-1 == bin_view) % ceiling
        [row,col] = find(ceiling == bin_view);
        sample_data_indep(k,4) = ceiling_x(row,col);
        sample_data_indep(k,5) = ceiling_y(row,col);
        sample_data_indep(k,6) = ceiling_z(row,col);
    elseif any(3203:3203+1280-1 == bin_view) % walls
        [row,col] = find(walls == bin_view);
        sample_data_indep(k,4) = walls_x(row,col);
        sample_data_indep(k,5) = walls_y(row,col);
        sample_data_indep(k,6) = walls_z(row,col);
    elseif any(4483:4483+160-1 == bin_view) % P1_BR
        [row,col] = find(P1_BR == bin_view);
        sample_data_indep(k,4) = P1_x(row,col);
        sample_data_indep(k,5) = P1_y(row,col);
        sample_data_indep(k,6) = PX_z(row,col);
    elseif any(4643:4643+160-1 == bin_view) % P2_BL
        [row,col] = find(P2_BL == bin_view);
        sample_data_indep(k,4) = P2_x(row,col);
        sample_data_indep(k,5) = P2_y(row,col);
        sample_data_indep(k,6) = PX_z(row,col);
    elseif any(4803:4803+160-1 == bin_view)
        [row,col] = find(P3_TR == bin_view); % P3_TR
        sample_data_indep(k,4) = P3_x(row,col);
        sample_data_indep(k,5) = P3_y(row,col);
        sample_data_indep(k,6) = PX_z(row,col);
    elseif any(4963:4963+160-1 == bin_view) % P4_TL
        [row,col] = find(P4_TL == bin_view);
        sample_data_indep(k,4) = P4_x(row,col);
        sample_data_indep(k,5) = P4_y(row,col);
        sample_data_indep(k,6) = PX_z(row,col);
    else % cue or hint
        invalid_rows = [invalid_rows k];
    end
    
end

sample_data_indep(invalid_rows,:) = [];
sample_data_dep(invalid_rows,:) = [];

%options_unc = optimoptions('fminunc','MaxFunctionEvaluations',10000,'FiniteDifferenceType','central','FunctionTolerance',1e-9,'StepTolerance',1e-10, 'MaxIterations', 1200);
options_unc = optimoptions('fminunc','SpecifyObjectiveGradient',true,'FiniteDifferenceType','central','FunctionTolerance',1e-9,'StepTolerance',1e-10, 'MaxIterations', 1200);
%options_unc = optimoptions('fminunc','SpecifyObjectiveGradient',true,'FiniteDifferenceType','central','Algorithm','quasi-newton','FunctionTolerance',1e-9,'StepTolerance',1e-10, 'MaxIterations', 1200);
options_con = optimoptions('fmincon','MaxFunctionEvaluations',10000,'SpecifyObjectiveGradient',true,'FiniteDifferenceType','central','FunctionTolerance',1e-9,'StepTolerance',1e-12, 'MaxIterations', 1200);
%options_con = optimoptions('fmincon','MaxFunctionEvaluations',10000,'FiniteDifferenceType','central','FunctionTolerance',1e-9,'StepTolerance',1e-10);

pillar_ub = 0.30 * tbin_size; % max allowed rate inside pillars, for constrained case

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Place model fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

place_params = cell(max_pfields, 1);
place_llh = cell(max_pfields, 1);

for pf = 1:max_pfields
    init_params = rand(pf*6,1) - 0.5;
    if nargin > 3 && pillcon
        A = []; B = []; Aeq = []; Beq = []; lb = []; ub = [];
        [temp_params, temp_llh] = fmincon(@(x)compute_llh_neg(x, pf, 0, sample_data_indep, sample_data_dep, 'place'), init_params, A, B, Aeq, Beq, lb, ub, @(x)pillar_cineq(x, pf, 0, pillar_a_x, pillar_a_y, pillar_ub, 'place'), options_con);
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
    init_params = 20*(rand(svf*10,1) - 0.5);
    [temp_params, temp_llh] = fminunc(@(x)compute_llh_neg(x, 0, svf, sample_data_indep, sample_data_dep, 'view'), init_params, options_unc);
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
        init_params = rand(pf*6+svf*10,1) - 0.5;
        [temp_params, temp_llh] = fminunc(@(x)compute_llh_neg(x, pf, svf, sample_data_indep, sample_data_dep, 'both'), init_params, options_unc);
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
eb_results.ThresVel = glm_data.ThresVel;
eb_results.UseMinObs = glm_data.UseMinObs;

eb_results.AIC_place = AIC_place;
eb_results.AIC_view = AIC_view;
eb_results.AIC_joint = AIC_joint;

if nargin > 3
    eb_results.pillcon = pillcon;
end


save('glm_eb_results_fields_1_5_1ms_unc.mat','eb_results','-v7.3');

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
    case 'view' % 10 parameters
        for k = 1:size(llh_entries,1)
            view_x = (data_indep(k,4) - 20) / 60;
            view_y = (data_indep(k,5) - 20) / 60;
            view_z = (data_indep(k,6) - 4) / 60;
            
            poiss_lambda = nan(svfield_num,1);
            for svf = 1:svfield_num
                exponent = input_params((svf-1)*10 + 1)*(view_x^2) + input_params((svf-1)*10 + 2)*(view_y^2) + input_params((svf-1)*10 + 3)*(view_z^2) + ...
                input_params((svf-1)*10 + 4)*view_x*view_y + input_params((svf-1)*10 + 5)*view_x*view_z + input_params((svf-1)*10 + 6)*view_y*view_z + ...
                input_params((svf-1)*10 + 7)*view_x + input_params((svf-1)*10 + 8)*view_y + input_params((svf-1)*10 + 9)*view_z + input_params((svf-1)*10 + 10);
                if exp(exponent) == 0
                    poiss_lambda(svf) = 1e-14;
                else
                    poiss_lambda(svf) = exp(exponent);
                end
                if nargout > 1 % gradient demanded
                    gradient(k, (svf-1)*10 + 1) = poiss_lambda(svf) * view_x^2;
                    gradient(k, (svf-1)*10 + 2) = poiss_lambda(svf) * view_y^2;
                    gradient(k, (svf-1)*10 + 3) = poiss_lambda(svf) * view_z^2;
                    gradient(k, (svf-1)*10 + 4) = poiss_lambda(svf) * view_x*view_y;
                    gradient(k, (svf-1)*10 + 5) = poiss_lambda(svf) * view_x*view_z;
                    gradient(k, (svf-1)*10 + 6) = poiss_lambda(svf) * view_y*view_z;
                    gradient(k, (svf-1)*10 + 7) = poiss_lambda(svf) * view_x;
                    gradient(k, (svf-1)*10 + 8) = poiss_lambda(svf) * view_y;
                    gradient(k, (svf-1)*10 + 9) = poiss_lambda(svf) * view_z;
                    gradient(k, (svf-1)*10 + 10) = poiss_lambda(svf);
                end
            end
            poiss_lambda = sum(poiss_lambda);
            if nargout > 1
                gradient(k, :) = ((data_dep(k)/poiss_lambda) - 1) * gradient(k, :);
            end
            llh_entries(k) = data_dep(k)*log(poiss_lambda) - poiss_lambda - log(factorial(data_dep(k)));
        end
    case 'both' % 6+10 parameters
        for k = 1:size(llh_entries,1)
            place_x = (data_indep(k,2) - 20) / 60;
            place_y = (data_indep(k,3) - 20) / 60;
            view_x = (data_indep(k,4) - 20) / 60;
            view_y = (data_indep(k,5) - 20) / 60;
            view_z = (data_indep(k,6) - 4) / 60;
            
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
                    view_params((svf-1)*10 + 1)*(view_x^2) + view_params((svf-1)*10 + 2)*(view_y^2) + view_params((svf-1)*10 + 3)*(view_z^2) + ...
                    view_params((svf-1)*10 + 4)*view_x*view_y + view_params((svf-1)*10 + 5)*view_x*view_z + view_params((svf-1)*10 + 6)*view_y*view_z + ...
                    view_params((svf-1)*10 + 7)*view_x + view_params((svf-1)*10 + 8)*view_y + view_params((svf-1)*10 + 9)*view_z + view_params((svf-1)*10 + 10);
                if exp(exponent) == 0
                    poiss_lambda(pfield_num + svf) = 1e-14;
                else
                    poiss_lambda(pfield_num + svf) = exp(exponent);
                end
                if nargout > 1
                    gradient(k, 6*pfield_num + (svf-1)*10 + 1) = poiss_lambda(pfield_num + svf) * view_x^2;
                    gradient(k, 6*pfield_num + (svf-1)*10 + 2) = poiss_lambda(pfield_num + svf) * view_y^2;
                    gradient(k, 6*pfield_num + (svf-1)*10 + 3) = poiss_lambda(pfield_num + svf) * view_z^2;
                    gradient(k, 6*pfield_num + (svf-1)*10 + 4) = poiss_lambda(pfield_num + svf) * view_x*view_y;
                    gradient(k, 6*pfield_num + (svf-1)*10 + 5) = poiss_lambda(pfield_num + svf) * view_x*view_z;
                    gradient(k, 6*pfield_num + (svf-1)*10 + 6) = poiss_lambda(pfield_num + svf) * view_y*view_z;
                    gradient(k, 6*pfield_num + (svf-1)*10 + 7) = poiss_lambda(pfield_num + svf) * view_x;
                    gradient(k, 6*pfield_num + (svf-1)*10 + 8) = poiss_lambda(pfield_num + svf) * view_y;
                    gradient(k, 6*pfield_num + (svf-1)*10 + 9) = poiss_lambda(pfield_num + svf) * view_z;
                    gradient(k, 6*pfield_num + (svf-1)*10 + 10) = poiss_lambda(pfield_num + svf);
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

place_params = input_params(1:6*pfield_num);
poiss_lambda = nan(pfield_num,length(x));
for pf = 1:pfield_num
    exponent = place_params((pf-1)*6 + 1)*scaled_x + place_params((pf-1)*6 + 2)*(scaled_x.^2) + ...
        place_params((pf-1)*6 + 3)*scaled_y + place_params((pf-1)*6 + 4)*(scaled_y.^2) + ...
        place_params((pf-1)*6 + 5)*scaled_x.*scaled_y + place_params((pf-1)*6 + 6);
    if exp(exponent) == 0
        poiss_lambda(pf,:) = 1e-14;
    else
        poiss_lambda(pf,:) = exp(exponent);
    end
end
poiss_lambda = sum(poiss_lambda(:,:),1);


c = poiss_lambda - upper_bound;
ceq = [];

end