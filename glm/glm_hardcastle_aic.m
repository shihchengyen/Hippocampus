function hc_results = glm_hardcastle_aic(glm_data, fc, smooth_params)
%
%   Reference: Hardcastle et al., 2017. A Multiplexed, Heterogeneous, and Adaptive Code for Navigation in Medial Entorhinal Cortex
%   Some code adapted from github.com/GiocomoLab/ln-model-of-mec-neurons
%
%   Be in directory of spiketrain.mat
%   Afterwards, do signrank tests, etc. with hardcastle_testing function
%   
%   PARAMETERS:
%   glm_data - uses either glm_vmpvData or glm_genData.
%   e.g. glm_hardcastle(glm_vmpvData(0.020), 10)
%
%   fc - number of folds used for cross-validation
%
%   smooth_params - 1x3 array of beta values for smoothing parameters,
%   in the order of [ place, headdirection, view ]. Optional argument,
%   defaults to [ 3, 3, 3 ] if not provided.

bin_stc = glm_data.bin_stc;
tbin_size = glm_data.tbin_size;
place_good_bins = glm_data.place_good_bins;
view_good_bins = glm_data.view_good_bins;

% Define bin geometry of the environment, in terms of:
% 1. no of floor_width bins, 2. no of wall_height bins, 3. no of
% pillar_height bins, 4. no of hd bins
floor_width = 40; wall_height = 8; pillar_height = 5; hd_bins = 60;
bin_geom = [ floor_width, wall_height, pillar_height, hd_bins ];
% additional intermediate variables for viewspace bins
viewbin_offset = 2; wall_perim = floor_width*4; pillar_width = floor_width/5; pillar_perim = pillar_width*4;

first_feature_bins = floor_width^2; % default 1600
second_feature_bins = hd_bins; % default 60
third_feature_bins = viewbin_offset + 2*floor_width^2 + wall_height*wall_perim + 4*pillar_height*pillar_perim; % default 5122

% Beta values for smoothing parameters
if exist('smooth_params', 'var')
    betas = smooth_params;
else
    betas = [ 3e0, 3e0, 3e0 ]; % [ beta_place, beta_hd, beta_view ], default value of 3 for each param
end

% Start to fill in x and y matrices
samples_total = size(bin_stc,1);
x = zeros(samples_total,first_feature_bins+second_feature_bins+third_feature_bins);
y = bin_stc(1:end,5);

for k = 1:samples_total
    x(k,bin_stc(k,2)) = 1;
    x(k,first_feature_bins+bin_stc(k,3)) = 1;
    x(k,first_feature_bins+second_feature_bins+bin_stc(k,4)) = 1;
end

% Filters for unoccupied place and view bins
place_filter = ones(first_feature_bins,1);
view_filter = ones(third_feature_bins,1);
place_filter(place_good_bins) = 0;
view_filter(view_good_bins) = 0;

%%% inputs done %%%

%%%%%%%%%% looping train-test splits %%%%%%%%%%

modelType = [1 1 1;
    1 1 0;
    1 0 1;
    0 1 1;
    1 0 0;
    0 1 0;
    0 0 1]; % testing for mixed model, ph model, pv model, hv model, place model, hd model, view model
modelName = {'phv', 'ph', 'pv', 'hv', 'place', 'headdirection', 'spatialview'};

folds = fc;
num_models = size(modelType, 1);

edges = round(linspace(1,length(y)+1, (5*folds)+1)); % splits dataset into 5xfold sections, to sample folds across time

testFit = nan(folds, num_models);
trainFit = nan(folds, num_models);
testFit_pure = nan(folds, num_models);
trainFit_pure = nan(folds, num_models);
paramsAll = cell(folds, num_models);

for model_type = 1:num_models % test different models on this dataset

    % Set bins that have no occurences to a value of -1e1 (sufficiently 
    % large negative number) during initialization of params
    bin_filter = [];
    if (modelType(model_type,1))
        bin_filter = [bin_filter; place_filter];
    end
    if (modelType(model_type,2))
        bin_filter = [bin_filter; ones(second_feature_bins,1)];
    end
    if (modelType(model_type,3))
        bin_filter = [bin_filter; view_filter];
    end
    bin_filter = find(bin_filter);
    
    % Random initialization of params for the first fold, then reuse
    % optimized params from the previous fold for subsequent folds
    param = 1e-3*randn(first_feature_bins*modelType(model_type,1) + second_feature_bins*modelType(model_type,2) + third_feature_bins*modelType(model_type,3), 1); % random initialization
    param(bin_filter) = -1e1; % set all bins that have no observations to -1e1 (sufficiently large negative number)
    
    disp(['Fitting model ', modelName{model_type}])

    for k = 1:folds
        fprintf('Fold %d of %d \n', k, folds)
%       disp('selectivity, params used, fold');
%       disp([cell_type model_type k]);

        test_ind  = [edges(k):edges(k+1)-1 edges(k+folds):edges(k+folds+1)-1 ...
            edges(k+2*folds):edges(k+2*folds+1)-1 edges(k+3*folds):edges(k+3*folds+1)-1 ...
            edges(k+4*folds):edges(k+4*folds+1)-1]; % grab indices for this fold, basically 5 smaller, spaced out sections

        train_ind = setdiff(1:length(y), test_ind); % remaining datapoints

        train_spikes = y(train_ind);
        test_spikes = y(test_ind);

        switch model_type
            case 1 % keep all columns (place, head direction and view info)
                train_A = x(train_ind,:);
                test_A = x(test_ind,:);
            case 2 % drop view info only
                train_A = x(train_ind,1:first_feature_bins+second_feature_bins);
                test_A = x(test_ind,1:first_feature_bins+second_feature_bins);
            case 3 % drop head direction info only
                train_A = x(train_ind,[1:first_feature_bins 1+first_feature_bins+second_feature_bins:end]);
                test_A = x(test_ind,[1:first_feature_bins 1+first_feature_bins+second_feature_bins:end]);
            case 4 % drop place info only
                train_A = x(train_ind,1+first_feature_bins:end);
                test_A = x(test_ind,1+first_feature_bins:end);
            case 5 % keep place info only
                train_A = x(train_ind,1:first_feature_bins);
                test_A = x(test_ind,1:first_feature_bins);
            case 6 % keep head direction info only
                train_A = x(train_ind,1+first_feature_bins:first_feature_bins+second_feature_bins);
                test_A = x(test_ind,1+first_feature_bins:first_feature_bins+second_feature_bins);
            case 7 % keep view info only
                train_A = x(train_ind,1+first_feature_bins+second_feature_bins:end);
                test_A = x(test_ind,1+first_feature_bins+second_feature_bins:end);            
        end

        opts = optimset('Gradobj','on','Hessian','on','Display','off');
        data{1} = train_A; 
        data{2} = train_spikes;
        init_param = param;
        % bottom part all adapted from reference github code
        [param] = fminunc(@(param) ln_poisson_model_vmpv(param,data,modelType(model_type,:),bin_geom,betas), init_param, opts);

        % test fit

        r = exp(test_A * param); 
        n = test_spikes; 
        meanFR_test = nanmean(test_spikes); 

        log_llh_test_model = nansum(r-n.*log(r)+log(factorial(n)))/sum(n); %note: log(gamma(n+1)) will be unstable if n is large (which it isn't here)
        log_llh_test_mean = nansum(meanFR_test-n.*log(meanFR_test)+log(factorial(n)))/sum(n);
        log_llh_test = (-log_llh_test_model + log_llh_test_mean);
        log_llh_test = log(2)*log_llh_test;
        
        aic_test = 2*sum(modelType(model_type))-2*log_llh_test;
        aic_test_model = 2*sum(modelType(model_type))-2*log(2)*(-log_llh_test_model);

        testFit(k, model_type) = aic_test;
        testFit_pure(k, model_type) = aic_test_model;

        % train fit

        r_train = exp(train_A * param); 
        n_train = train_spikes; 
        meanFR_train = nanmean(train_spikes);   

        log_llh_train_model = nansum(r_train-n_train.*log(r_train)+log(factorial(n_train)))/sum(n_train);
        log_llh_train_mean = nansum(meanFR_train-n_train.*log(meanFR_train)+log(factorial(n_train)))/sum(n_train);
        log_llh_train = (-log_llh_train_model + log_llh_train_mean);
        log_llh_train = log(2)*log_llh_train;

        aic_train = 2*sum(modelType(model_type))-2*log_llh_train;
        aic_train_model = 2*sum(modelType(model_type))-2*log(2)*(-log_llh_train_model);
        
        trainFit(k, model_type) = aic_train;
        trainFit_pure(k, model_type) = aic_train_model;

        paramsAll{k, model_type} = param;

    end
    
end
    
%glm_hardcastle_results.inputs = [x y];
hc_results.training_fits = trainFit; % for each dataset (phv, ph, pv, hv, place, hd, view behavior), store n_folds x model_type training likelihood values (with subtraction of mean model)
hc_results.training_fits_pure = trainFit_pure; % training likelihood without comparison to mean model
hc_results.testing_fits = testFit; % test likelihood with comparison to mean model
hc_results.testing_fits_pure = testFit_pure; % test likelihood without comparison to mean model
hc_results.params_consol = paramsAll; % weights stored here

hc_results.tbin_size = tbin_size;
hc_results.num_folds = fc;
hc_results.smoothing_beta = betas;
%hc_results.ThresVel = glm_data.ThresVel;
%hc_results.UseMinObs = glm_data.UseMinObs;

save('glm_hardcastle_results.mat','hc_results','-v7.3');

end

