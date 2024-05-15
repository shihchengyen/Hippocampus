function hardcastle_testing(data, use_training, use_pure)

    if ~exist('use_training', 'var')
        use_training = 0;
    end
    if ~exist('use_demeaned', 'var')
        use_pure = 0;
    end
    
    if use_training == 1
        mean_fits = data.training_fits;
        if use_pure == 1
            fits = data.training_fits_pure;
        else
            fits = mean_fits;
        end
    else
        mean_fits = data.testing_fits;
        if use_pure == 1
            fits = data.testing_fits_pure;
        else
            fits = mean_fits;
        end
    end

    % 1: PHV, 2: PH, 3: PV, 4: HV, 5: P, 6: H, 7: V
    
    disp('Single variable model tests');
    disp(' ');
    
    x = signrank(fits(:,5), fits(:,6), 'tail', 'right');
    disp(['P > H: ' num2str(x)]);
    x = signrank(fits(:,6), fits(:,5), 'tail', 'right');
    disp(['H > P: ' num2str(x)]);  
    disp(' ');
    
    x = signrank(fits(:,5), fits(:,7), 'tail', 'right');
    disp(['P > V: ' num2str(x)]);
    x = signrank(fits(:,7), fits(:,5), 'tail', 'right');
    disp(['V > P: ' num2str(x)]);   
    disp(' ');
    
    x = signrank(fits(:,6), fits(:,7), 'tail', 'right');
    disp(['H > V: ' num2str(x)]);
    x = signrank(fits(:,7), fits(:,6), 'tail', 'right');
    disp(['V > H: ' num2str(x)]);   
    disp(' ');
    
    disp('Two variable model tests');
    disp(' ');
    
    x = signrank(fits(:,2), fits(:,5), 'tail', 'right');
    disp(['PH > P: ' num2str(x)]);
    x = signrank(fits(:,2), fits(:,6), 'tail', 'right');
    disp(['PH > H: ' num2str(x)]);
    disp(' ');
    x = signrank(fits(:,3), fits(:,5), 'tail', 'right');
    disp(['PV > P: ' num2str(x)]);
    x = signrank(fits(:,3), fits(:,7), 'tail', 'right');
    disp(['PV > V: ' num2str(x)]);
    disp(' ');
    x = signrank(fits(:,4), fits(:,6), 'tail', 'right');
    disp(['HV > H: ' num2str(x)]);
    x = signrank(fits(:,4), fits(:,7), 'tail', 'right');
    disp(['HV > V: ' num2str(x)]);
    disp(' ');
    
    disp('Three variable model tests');
    disp(' ');
    
    x = signrank(fits(:,1), fits(:,2), 'tail', 'right');
    disp(['PHV > PH: ' num2str(x)]);
    x = signrank(fits(:,1), fits(:,3), 'tail', 'right');
    disp(['PHV > PV: ' num2str(x)]);
    x = signrank(fits(:,1), fits(:,4), 'tail', 'right');
    disp(['PHV > HV: ' num2str(x)]);
    disp(' ');
    
    disp('Mean rate tests');
    disp(' ');
    
    x = signrank(mean_fits(:,1), 0, 'tail', 'right');
    disp(['PHV mean rate test: ' num2str(x)]);
    x = signrank(mean_fits(:,2), 0, 'tail', 'right');
    disp([' PH mean rate test: ' num2str(x)]);
    x = signrank(mean_fits(:,3), 0, 'tail', 'right');
    disp([' PV mean rate test: ' num2str(x)]);
    x = signrank(mean_fits(:,4), 0, 'tail', 'right');
    disp([' HV mean rate test: ' num2str(x)]);
    x = signrank(mean_fits(:,5), 0, 'tail', 'right');
    disp(['  P mean rate test: ' num2str(x)]);
    x = signrank(mean_fits(:,6), 0, 'tail', 'right');
    disp(['  H mean rate test: ' num2str(x)]);
    x = signrank(mean_fits(:,7), 0, 'tail', 'right');
    disp(['  V mean rate test: ' num2str(x)]);
    
end