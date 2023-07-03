function hc_results = glm_hardcastle(glm_data, fc)
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
%   fc - refers to the number of folds used for cross-validation
%


bin_stc = glm_data.bin_stc;
tbin_size = glm_data.tbin_size;

first_feature_bins = 1600;
second_feature_bins = 5122;

% Beta values for smoothing parameters
beta_place = 3e0;
beta_view = 3e0;

% Start to fill in x and y matrices
samples_total = size(bin_stc,1);
x = zeros(samples_total,first_feature_bins+second_feature_bins);
y = zeros(samples_total,1);

for k = 1:samples_total
    x(k,bin_stc(k,2)) = 1;
    x(k,first_feature_bins+bin_stc(k,3)) = 1;
end

spike_latest = 1;
for k = 1:size(bin_stc,1)
    for g = spike_latest:length(spiketimes)
        if spiketimes(g) > bin_stc(k,1) + tbin_size
            spike_latest = g;
            break
        elseif spiketimes(g) > bin_stc(k,1)
            y(k) = y(k) + 1;
        end
    end
end

%%% inputs done %%%

%%%%%%%%%% looping train-test splits %%%%%%%%%%

modelType = [1 1; 
    1 0; 
    0 1]; % testing for mixed model, place model, view model

folds = fc;

edges = round(linspace(1,length(y)+1, (5*folds)+1)); % splits dataset into 5xfold sections, to sample folds across time

testFit = nan(folds, 3);
trainFit = nan(folds, 3);
testFit_pure = nan(folds, 3);
trainFit_pure = nan(folds, 3);
paramsAll = cell(folds, 3);

for model_type = 1:3 % test different models on this dataset

    param = 1e-3*randn(first_feature_bins*modelType(model_type,1) + second_feature_bins*modelType(model_type,2), 1); % random initialization
    disp(['Fitting model ' num2str(model_type)])

    for k = 1:folds
%                 disp('selectivity, params used, fold');
        param = 1e-3*randn(first_feature_bins*modelType(model_type,1) + second_feature_bins*modelType(model_type,2), 1); % reset across folds
%                 disp([cell_type model_type k]);

        disp(num2str(k))

        test_ind  = [edges(k):edges(k+1)-1 edges(k+folds):edges(k+folds+1)-1 ...
            edges(k+2*folds):edges(k+2*folds+1)-1 edges(k+3*folds):edges(k+3*folds+1)-1 ...
            edges(k+4*folds):edges(k+4*folds+1)-1]; % grab indices for this fold, basically 5 smaller, spaced out sections

        train_ind = setdiff(1:length(y), test_ind); % remaining datapoints

        train_spikes = y(train_ind);
        test_spikes = y(test_ind);

        if model_type == 1 % keep all columns (both place and view info)
            train_A = x(train_ind,:);
            test_A = x(test_ind,:);
        elseif model_type == 2 % keep place info only
            train_A = x(train_ind,1:first_feature_bins);
            test_A = x(test_ind,1:first_feature_bins);
        else % keep view info only
            train_A = x(train_ind,1+first_feature_bins:end);
            test_A = x(test_ind,1+first_feature_bins:end);            
        end

        opts = optimset('Gradobj','on','Hessian','on','Display','off');
        data{1} = train_A; 
        data{2} = train_spikes;
        init_param = param;
        % bottom part all adapted from reference github code
        [param] = fminunc(@(param) ln_poisson_model_vmpv(param,data,modelType(model_type,:),first_feature_bins,second_feature_bins,beta_place,beta_view), init_param, opts);

        % test fit

        r = exp(test_A * param); 
        n = test_spikes; 
        meanFR_test = nanmean(test_spikes); 

        log_llh_test_model = nansum(r-n.*log(r)+log(factorial(n)))/sum(n); %note: log(gamma(n+1)) will be unstable if n is large (which it isn't here)
        log_llh_test_mean = nansum(meanFR_test-n.*log(meanFR_test)+log(factorial(n)))/sum(n);
        log_llh_test = (-log_llh_test_model + log_llh_test_mean);
        log_llh_test = log(2)*log_llh_test;

        testFit(k, model_type) = log_llh_test;
        testFit_pure(k, model_type) = log(2)*(-log_llh_test_model);

        % train fit

        r_train = exp(train_A * param); 
        n_train = train_spikes; 
        meanFR_train = nanmean(train_spikes);   

        log_llh_train_model = nansum(r_train-n_train.*log(r_train)+log(factorial(n_train)))/sum(n_train);
        log_llh_train_mean = nansum(meanFR_train-n_train.*log(meanFR_train)+log(factorial(n_train)))/sum(n_train);
        log_llh_train = (-log_llh_train_model + log_llh_train_mean);
        log_llh_train = log(2)*log_llh_train;

        trainFit(k, model_type) = log_llh_train;
        trainFit_pure(k, model_type) = log(2)*(-log_llh_train_model);

        paramsAll{k, model_type} = param;

    end
    
end
    
%glm_hardcastle_results.inputs = [x y];
hc_results.training_fits = trainFit; % for each dataset (joint, place, view behavior), store n_folds x model_type training likelihood values (with subtraction of mean model)
hc_results.training_fits_pure = trainFit_pure; % training likelihood without comparison to mean model
hc_results.testing_fits = testFit; % test likelihood with comparison to mean model
hc_results.testing_fits_pure = testFit_pure; % test likelihood without comparison to mean model
hc_results.params_consol = paramsAll; % weights stored here

hc_results.smoothing_beta = cell(2,1);
hc_results.smoothing_beta{1} = beta_place;
hc_results.smoothing_beta{2} = beta_view;

hc_results.tbin_size = tbin_size;
%hc_results.ThresVel = glm_vmpvData.ThresVel;
%hc_results.UseMinObs = glm_vmpvData.UseMinObs;

save(['glm_hardcastle_results_' num2str(tbin_size) '.mat'],'hc_results','-v7.3');

end

