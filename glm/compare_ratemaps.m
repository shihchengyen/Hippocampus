function [similarity_scores, hc_results] = compare_ratemaps(hc_results, metric, reference)
    % Compares model-fitted ratemap to actual ratemap from data. Requires
    % vmpc/vmhd/vmsv objects to be in the current directory before running.
    %
    % PARAMETERS:
    % hc_results - struct output of glm_hardcastle, after running
    % classification via select_best_model
    % metric - function to be used in comparing ratemaps OR string of an
    % already implemented function, should take in 3 arguments (arr1:
    % numeric array, arr2: numeric array, optArgs: cell array of optional
    % args) and output a single scalar value
    %
    % Current format of optArgs: optArgs{1} - dur_map (for use in weighted
    % cosine similarity)
    
    % List of implemented functions:
    implemented_funcs = {'cosine_similarity', 'weighted_cosine_similarity'};

    params = hc_results.params_consol;
    tbin_size = hc_results.tbin_size;
    num_folds = hc_results.num_folds;
    
    global num_place; global num_hd; global num_view;
    num_place = 1600;
    num_hd = 60;
    num_view = 5122;
    
    if ischar(metric)
        if ismember(metric, implemented_funcs)
            metric = str2func(metric);
        else
            error(['Invalid function name provided: ' metric])
        end
    end
    
    if ~exist('reference', 'var')
        reference = 'smoothed';  % use smoothed ratemaps from vmobj by default
    end        

    % load in ratemaps from actual data
    switch reference
        case 'raw'
            % raw ratemaps from vmobj
            load('vmpc.mat', 'vmp');
            load('vmhd.mat', 'vmd');
            load('vmsv.mat', 'vms');
            
            actual_place = repmat({vmp.data.maps_raw'}, num_folds, 1);
            actual_hd = repmat({vmd.data.maps_raw'}, num_folds, 1);
            actual_view = repmat({vms.data.maps_raw'}, num_folds, 1);
            
            place_dur = repmat({vmp.data.dur_raw'}, num_folds, 1);
            hd_dur = repmat({vmd.data.dur_raw'}, num_folds, 1);
            view_dur = repmat({vms.data.dur_raw'}, num_folds, 1);
            
        case 'smoothed'
            % smoothed ratemaps from vmobj
            load('vmpc.mat', 'vmp');
            load('vmhd.mat', 'vmd');
            load('vmsv.mat', 'vms');
            
            actual_place = repmat({vmp.data.maps_adsm'}, num_folds, 1);
            actual_hd = repmat({vmd.data.maps_sm'}, num_folds, 1);
            actual_view = repmat({vms.data.maps_adsm'}, num_folds, 1);
            
            place_dur = repmat({vmp.data.dur_raw'}, num_folds, 1);
            hd_dur = repmat({vmd.data.dur_raw'}, num_folds, 1);
            view_dur = repmat({vms.data.dur_raw'}, num_folds, 1);
            
        case 'training'
            % raw ratemaps from training fold
            if exist('vmpvData.mat', 'file')
                load('vmpvData.mat', 'vmpvData');
            else
                vmpvData = glm_vmpvData(tbin_size);
            end
            bin_stc = vmpvData.bin_stc;
            
            actual_place = repmat({zeros(num_place, 1)}, num_folds, 1);
            actual_hd = repmat({zeros(num_hd, 1)}, num_folds, 1);
            actual_view = repmat({zeros(num_view, 1)}, num_folds, 1);
            
            place_dur = repmat({zeros(num_place, 1)}, num_folds, 1);
            hd_dur = repmat({zeros(num_hd, 1)}, num_folds, 1);
            view_dur = repmat({zeros(num_view, 1)}, num_folds, 1);
            
            fold_edges = round(linspace(1,size(bin_stc, 1)+1, (5*num_folds)+1));
            for k = 1:num_folds
                % Get slice indices for training data for the current fold
                test_ind  = [fold_edges(k):fold_edges(k+1)-1 fold_edges(k+num_folds):fold_edges(k+num_folds+1)-1 ...
                    fold_edges(k+2*num_folds):fold_edges(k+2*num_folds+1)-1 fold_edges(k+3*num_folds):fold_edges(k+3*num_folds+1)-1 ...
                    fold_edges(k+4*num_folds):fold_edges(k+4*num_folds+1)-1];
                train_ind = setdiff(1:size(bin_stc, 1), test_ind);
                bin_stc_fold = bin_stc(train_ind,:);
                
                for i = 1:size(bin_stc_fold, 1)
                    % Fill in spike count and duration maps for the current fold
                    actual_place{k}(bin_stc_fold(i,2)) = actual_place{k}(bin_stc_fold(i,2)) + bin_stc_fold(i,5);
                    actual_hd{k}(bin_stc_fold(i,3)) = actual_hd{k}(bin_stc_fold(i,3)) + bin_stc_fold(i,5);
                    actual_view{k}(bin_stc_fold(i,4)) = actual_view{k}(bin_stc_fold(i,4)) + bin_stc_fold(i,5);

                    place_dur{k}(bin_stc_fold(i,2)) = place_dur{k}(bin_stc_fold(i,2)) + tbin_size;
                    hd_dur{k}(bin_stc_fold(i,3)) = hd_dur{k}(bin_stc_fold(i,3)) + tbin_size;
                    view_dur{k}(bin_stc_fold(i,4)) = view_dur{k}(bin_stc_fold(i,4)) + tbin_size;
                end
                
                % Calculate average firing rate per bin, replace divide-by-zero results with NaN
                actual_place{k} = actual_place{k} ./ place_dur{k};
                actual_place{k}(isinf(actual_place{k})) = NaN;
                actual_hd{k} = actual_hd{k} ./ hd_dur{k};
                actual_hd{k}(isinf(actual_hd{k})) = NaN;
                actual_view{k} = actual_view{k} ./ view_dur{k};
                actual_view{k}(isinf(actual_view{k})) = NaN;
                
            end
    end
    
    % create 1x7 cell array to store all results
    similarity_scores = cell(7,1);
    
    for model = 1:7
        % extract relevant parameters for each variable space
        switch model
            case 1  % phv
                model_params = cell2mat(params(:, 1)')';
                place_params = model_params(:, 1:num_place);
                hd_params = model_params(:, num_place+1:num_place+num_hd);
                view_params = model_params(:, num_place+num_hd+1:end);

            case 2  % ph
                model_params = cell2mat(params(:, 2)')';
                place_params = model_params(:, 1:num_place);
                hd_params = model_params(:, num_place+1:end);

            case 3  % pv
                model_params = cell2mat(params(:, 3)')';
                place_params = model_params(:, 1:num_place);
                view_params = model_params(:, num_place+1:end);

            case 4  % hv
                model_params = cell2mat(params(:, 4)')';
                hd_params = model_params(:, 1:num_hd);
                view_params = model_params(:, num_hd+1:end);

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
                place_similarity(fc) = metric(model_place, actual_place{fc}, {place_dur{fc}});
            end
        end

        if model == 1 || model == 2 || model == 4 || model == 6  % headdirection/ph/hv/phv
            hd_similarity = nan(num_folds, 1);
            for fc = 1:num_folds
                model_hd = nan(60,1);
                for k = 1:size(model_hd,1)
                    model_hd(k) = exp(hd_params(fc, k))/tbin_size;
                end
                hd_similarity(fc) = metric(model_hd, actual_hd{fc}, {hd_dur{fc}});
            end
        end

        if model == 1 || model == 3 || model == 4 || model == 7  % spatialview/pv/hv/phv
            view_similarity = nan(num_folds, 1);
            for fc = 1:num_folds
                model_view = nan(5122,1);
                for k = 1:size(model_view,1)
                    model_view(k) = exp(view_params(fc, k))/tbin_size;
                end
                view_similarity(fc) = metric(model_view, actual_view{fc}, {view_dur{fc}});
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

function score = cosine_similarity(arr1, arr2, optArgs)
    score = nansum(arr1.*arr2)/(sqrt(nansum(arr1.^2))*sqrt(nansum(arr2.^2)));
end

function score = weighted_cosine_similarity(arr1, arr2, optArgs)
    dur_map = optArgs{1};
    weight = dur_map / nansum(dur_map);
    score = nansum(weight.*arr1.*arr2)/(sqrt(nansum(arr1.^2))*sqrt(nansum(arr2.^2)));
end
