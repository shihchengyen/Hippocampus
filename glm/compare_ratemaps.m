function [similarity_scores, hc_results] = compare_ratemaps(hc_results, metric)
    % Compares model-fitted ratemap to actual ratemap from data. Requires
    % vmpc/vmhd/vmsv objects to be in the current directory before running.
    %
    % PARAMETERS:
    % hc_results - struct output of glm_hardcastle, after running
    % classification via select_best_model
    % metric - function to be used in comparing ratemaps OR string of an
    % already implemented function, should output a single scalar value
    
    % List of implemented functions:
    implemented_funcs = {'cosine_similarity', 'weighted_cosine_similarity'};

    params = hc_results.params_consol;
    tbin_size = hc_results.tbin_size;
    num_folds = hc_results.num_folds;
    
    if ischar(metric) && ismember(metric, implemented_funcs)
        metric = str2func(metric);
    end

    % load in ratemaps from actual data
    global vmp; global vmd; global vms;
    load('vmpc.mat', 'vmp');
    load('vmhd.mat', 'vmd');
    load('vmsv.mat', 'vms');
    actual_place = vmp.data.maps_adsm';
    actual_hd = vmd.data.maps_sm';
    actual_view = vms.data.maps_adsm';
    
    % create 1x7 cell array to store all results
    similarity_scores = cell(7,1);
    
    for model = 1:7
        % extract relevant parameters for each variable space
        switch model
            case 1  % phv
                model_params = cell2mat(params(:, 1)')';
                place_params = model_params(:, 1:1600);
                hd_params = model_params(:, 1601:1600+60);
                view_params = model_params(:, 1601+60:1600+60+5122);

            case 2  % ph
                model_params = cell2mat(params(:, 2)')';
                place_params = model_params(:, 1:1600);
                hd_params = model_params(:, 1601:1600+60);

            case 3  % pv
                model_params = cell2mat(params(:, 3)')';
                place_params = model_params(:, 1:1600);
                view_params = model_params(:, 1601:1600+5122);

            case 4  % hv
                model_params = cell2mat(params(:, 4)')';
                hd_params = model_params(:, 1:60);
                view_params = model_params(:, 61:60+5122);

            case 5  % place
                place_params = cell2mat(params(:, 5)')';
                
            case 6  % headdirection
                hd_params = cell2mat(params(:, 6)')';

            case 7  % spatialview
                view_params = cell2mat(params(:, 7)')';
        end

        % compute similarity scores between model-fitted and actual ratemaps
        % for each variable
        if model == 1 || model == 2 || model == 3 || model == 5  % place/ph/pv/phv
            place_similarity = nan(num_folds, 1);
            for fc = 1:num_folds
                model_place = nan(1600,1);
                for k = 1:size(model_place,1)
                    model_place(k) = exp(place_params(fc, k))/tbin_size;
                end
                place_similarity(fc) = metric(model_place, actual_place);
            end
        end

        if model == 1 || model == 2 || model == 4 || model == 6  % headdirection/ph/hv/phv
            hd_similarity = nan(num_folds, 1);
            for fc = 1:num_folds
                model_hd = nan(60,1);
                for k = 1:size(model_hd,1)
                    model_hd(k) = exp(hd_params(fc, k))/tbin_size;
                end
                hd_similarity(fc) = metric(model_hd, actual_hd);
            end
        end

        if model == 1 || model == 3 || model == 4 || model == 7  % spatialview/pv/hv/phv
            view_similarity = nan(num_folds, 1);
            for fc = 1:num_folds
                model_view = nan(5122,1);
                for k = 1:size(model_view,1)
                    model_view(k) = exp(view_params(fc, k))/tbin_size;
                end
                view_similarity(fc) = metric(model_view, actual_view);
            end
        end

        % save results back to hc_results as a new data field
        switch model
            case 1
                similarity_scores{model} = [place_similarity, hd_similarity, view_similarity];
            case 2
                similarity_scores{model} = [place_similarity, hd_similarity];
            case 3
                similarity_scores{model} = [place_similarity, view_similarity];
            case 4
                similarity_scores{model} = [hd_similarity, view_similarity];
            case 5
                similarity_scores{model} = place_similarity;
            case 6
                similarity_scores{model} = hd_similarity;
            case 7
                similarity_scores{model} = view_similarity;
        end
    end
    
    hc_results.similarity_scores = similarity_scores;
    save('glm_hardcastle_results.mat','hc_results','-v7.3');
    
end

function score = cosine_similarity(arr1, arr2)
    score = nansum(arr1.*arr2)/(sqrt(nansum(arr1.^2))*sqrt(nansum(arr2.^2)));
end

function score = weighted_cosine_similarity(arr1, arr2)
    switch length(arr1)
        case 1600  % place
            global vmp;
            dur_map = vmp.data.dur_raw;
        case 60  % hd
            global vmd;
            dur_map = vmd.data.dur_raw;
        case 5122  % view
            global vms;
            dur_map = vms.data.dur_raw;
    end
    weight = dur_map / nansum(dur_map);
    score = nansum(weight.*arr1.*arr2)/(sqrt(nansum(arr1.^2))*sqrt(nansum(arr2.^2)));
end
