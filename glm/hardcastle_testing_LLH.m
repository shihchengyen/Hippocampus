function hardcastle_testing_LLH(data, use_training, use_demeaned)

    if ~exist('use_training', 'var')
        use_training = 0;
    end
    if ~exist('use_demeaned', 'var')
        use_demeaned = 0;
    end
    
    if use_training == 1
        mean_fits = data.training_fits;
        if use_demeaned == 0
            fits = data.training_fits_pure;
        else
            fits = mean_fits;
        end
    else
        mean_fits = data.testing_fits;
        if use_demeaned == 0
            fits = data.testing_fits_pure;
        else
            fits = mean_fits;
        end
    end

    % 1: PHV, 2: PH, 3: PV, 4: HV, 5: P, 6: H, 7: V
    model_names = {'PHV', 'PH', 'PV', 'HV', 'P', 'H', 'V'};
    mean_LLH_values = nanmean(fits);
    
    disp('Single variable model LLHs');
    
    x = mean_LLH_values(5); y = signrank(mean_fits(:,5), 0, 'tail', 'right');
    disp(['P model LLH: ' num2str(x) ', mean rate test: ' num2str(y)]);
    x = mean_LLH_values(6); y = signrank(mean_fits(:,6), 0, 'tail', 'right');
    disp(['H model LLH: ' num2str(x) ', mean rate test: ' num2str(y)]);
    x = mean_LLH_values(7); y = signrank(mean_fits(:,7), 0, 'tail', 'right');
    disp(['V model LLH: ' num2str(x) ', mean rate test: ' num2str(y)]);
    
    singleModels = 5:7;
    [~,top1] = max(mean_LLH_values(singleModels)); top1 = top1 + singleModels(1)-1;
    disp(['Best single variable model: ' model_names{top1}]);
    disp(' ');
    
    disp('Double variable model LLHs');
    
    % find the best double model that includes the single model
    if top1 == 5 % P -> PH, PV
        doubleModels = [2 3];      
    elseif top1 == 6 % H -> PH, HV
        doubleModels = [2 4];
    else % V -> PV, HV
        doubleModels = [3 4];
    end
    
    for i = 1:length(doubleModels)
        m = doubleModels(i);
        x = mean_LLH_values(m); y = signrank(mean_fits(:,m), 0, 'tail', 'right');
        disp([model_names{m} ' model LLH: ' num2str(x) ', mean rate test: ' num2str(y)]);
    end
    
    [~,top2] = max(mean_LLH_values(doubleModels));
    top2 = doubleModels(top2);
    disp(['Best double variable model: ' model_names{top2}]);
    
    x = signrank(fits(:,top2), fits(:,top1), 'tail', 'right');
    disp([model_names{top2} ' > ' model_names{top1} ': ' num2str(x)]);
    disp(' ');
    
    disp('Triple variable model LLHs');
    
    top3 = 1;
    x = mean_LLH_values(top3); y = signrank(mean_fits(:,top3), 0, 'tail', 'right');
    disp([model_names{top3} ' model LLH: ' num2str(x) ', mean rate test: ' num2str(y)]);
    x = signrank(fits(:,top3), fits(:,top2), 'tail', 'right');
    disp([model_names{top3} ' > ' model_names{top2} ': ' num2str(x)]);
    disp(' ');
    
end