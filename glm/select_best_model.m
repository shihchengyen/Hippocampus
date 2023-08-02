%% Description
% Code adapted from: https://github.com/GiocomoLab/ln-model-of-mec-neurons/blob/master/select_best_model.m

% This function implements forward feature selection in order to determine
% the simplest model that best describes neural spiking. First, the
% highest-performing single-variable model is identified. Then, the
% highest-performing double-variable model that includes the
% single-variable model is identified. This continues until the full model
% is identified. Next, statistical tests are applied to see if including
% extra variables significantly improves model performance. The first time
% that including variable does NOT signficantly improve performance, the
% procedure is stopped and the model at that point is recorded as the
% selected model.

% the model indexing scheme:
% phv, ph, pv, hv, p,  h,  v
%  1   2   3   4   5   6   7

function selected_model = select_best_model(hc_results, p_sig)

    if ~exist('p_sig', 'var')
        p_sig = 0.05;
    end

    LLH_values = hc_results.testing_fits;

    % find the best single model
    singleModels = 5:7;
    [~,top1] = max(nanmean(LLH_values(:,singleModels))); top1 = top1 + singleModels(1)-1;

    % find the best double model that includes the single model
    if top1 == 5 % P -> PH, PV
        [~,top2] = max(nanmean(LLH_values(:,[2 3])));
        vec = [2 3]; top2 = vec(top2);
    elseif top1 == 6 % H -> PH, HV
        [~,top2] = max(nanmean(LLH_values(:,[2 4])));
        vec = [2 4]; top2 = vec(top2);
    else % V -> PV, HV
        [~,top2] = max(nanmean(LLH_values(:,[3 4])));
        vec = [3 4]; top2 = vec(top2);
    end

    top3 = 1;
    LLH1 = LLH_values(:,top1); LLH2 = LLH_values(:,top2); LLH3 = LLH_values(:,top3);

    [p_llh_12,~] = signrank(LLH2,LLH1,'tail','right');
    [p_llh_23,~] = signrank(LLH3,LLH2,'tail','right');

    if p_llh_12 < p_sig % double model is sig. better
        if p_llh_23 < p_sig  % full model is sig. better
            selected_model = top3; % full model
        else
            selected_model = top2; % double model
        end
    else
        selected_model = top1; % single model
    end

    % re-set if selected model is not above baseline
    pval_baseline = signrank(LLH_values(:,selected_model),[],'tail','right');

    if pval_baseline > p_sig
        selected_model = 8; % unclassified
    end
    
    model_names = {'phv', 'ph', 'pv', 'hv',
        'place', 'headdirection', 'spatialview', 'unclassified'};
    selected_model = model_names{selected_model};
    
end
