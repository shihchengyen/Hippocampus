function hardcastle_testing(data, use_training, use_demeaned)

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
    
    disp('Single variable model tests');
    disp(' ');
    
    x = signrank(fits(:,5), fits(:,6), 'tail', 'both');
    disp(['P vs H, 2 tail: ' num2str(x)]);
    x = signrank(fits(:,5), fits(:,6), 'tail', 'right');
    disp(['P > H, 1 tail: ' num2str(x)]);
    x = signrank(fits(:,6), fits(:,5), 'tail', 'right');
    disp(['H > P, 1 tail: ' num2str(x)]);  
    disp(' ');
    
    x = signrank(fits(:,5), fits(:,7), 'tail', 'both');
    disp(['P vs V, 2 tail: ' num2str(x)]);
    x = signrank(fits(:,5), fits(:,7), 'tail', 'right');
    disp(['P > V, 1 tail: ' num2str(x)]);
    x = signrank(fits(:,7), fits(:,5), 'tail', 'right');
    disp(['V > P, 1 tail: ' num2str(x)]);   
    disp(' ');
    
    x = signrank(fits(:,6), fits(:,7), 'tail', 'both');
    disp(['H vs V, 2 tail: ' num2str(x)]);
    x = signrank(fits(:,6), fits(:,7), 'tail', 'right');
    disp(['H > V, 1 tail: ' num2str(x)]);
    x = signrank(fits(:,7), fits(:,6), 'tail', 'right');
    disp(['V > H, 1 tail: ' num2str(x)]);   
    disp(' ');
    
    disp('Two variable model tests');
    disp(' ');
    
    x = signrank(fits(:,2), fits(:,5), 'tail', 'right');
    disp(['PH > P, 1 tail: ' num2str(x)]);
    x = signrank(fits(:,2), fits(:,6), 'tail', 'right');
    disp(['PH > H, 1 tail: ' num2str(x)]);
    disp(' ');
    x = signrank(fits(:,3), fits(:,5), 'tail', 'right');
    disp(['PV > P, 1 tail: ' num2str(x)]);
    x = signrank(fits(:,3), fits(:,7), 'tail', 'right');
    disp(['PV > V, 1 tail: ' num2str(x)]);
    disp(' ');
    x = signrank(fits(:,4), fits(:,6), 'tail', 'right');
    disp(['HV > H, 1 tail: ' num2str(x)]);
    x = signrank(fits(:,4), fits(:,7), 'tail', 'right');
    disp(['HV > V, 1 tail: ' num2str(x)]);
    disp(' ');
    
    disp('Three variable model tests');
    disp(' ');
    
    x = signrank(fits(:,1), fits(:,2), 'tail', 'right');
    disp(['PHV > PH, 1 tail: ' num2str(x)]);
    x = signrank(fits(:,1), fits(:,3), 'tail', 'right');
    disp(['PHV > PV, 1 tail: ' num2str(x)]);
    x = signrank(fits(:,1), fits(:,4), 'tail', 'right');
    disp(['PHV > HV, 1 tail: ' num2str(x)]);
    disp(' ');
    
    disp('Mean rate tests');
    disp(' ');
    
    x = signrank(mean_fits(:,1), 0, 'tail', 'right');
    disp(['PHV mean rate test, 1 tail: ' num2str(x)]);
    x = signrank(mean_fits(:,2), 0, 'tail', 'right');
    disp([' PH mean rate test, 1 tail: ' num2str(x)]);
    x = signrank(mean_fits(:,3), 0, 'tail', 'right');
    disp([' PV mean rate test, 1 tail: ' num2str(x)]);
    x = signrank(mean_fits(:,4), 0, 'tail', 'right');
    disp([' HV mean rate test, 1 tail: ' num2str(x)]);
    x = signrank(mean_fits(:,5), 0, 'tail', 'right');
    disp(['  P mean rate test, 1 tail: ' num2str(x)]);
    x = signrank(mean_fits(:,6), 0, 'tail', 'right');
    disp(['  H mean rate test, 1 tail: ' num2str(x)]);
    x = signrank(mean_fits(:,7), 0, 'tail', 'right');
    disp(['  V mean rate test, 1 tail: ' num2str(x)]);
    
end