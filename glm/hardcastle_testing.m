function hardcastle_testing(data, use_training, use_demeaned)

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

    x = signrank(fits(:,2), fits(:,3), 'tail', 'both');
    disp(['P vs V, 2 tail: ' num2str(x)]);
    x = signrank(fits(:,2), fits(:,3), 'tail', 'right');
    disp(['P > V, 1 tail: ' num2str(x)]);
    x = signrank(fits(:,3), fits(:,2), 'tail', 'right');
    disp(['V > P, 1 tail: ' num2str(x)]);   
    disp(' ');
    x = signrank(fits(:,1), fits(:,2), 'tail', 'right');
    disp(['J > P, 1 tail: ' num2str(x)]);
    x = signrank(fits(:,1), fits(:,3), 'tail', 'right');
    disp(['J > V, 1 tail: ' num2str(x)]);
    disp(' ');
    x = signrank(mean_fits(:,1), 0, 'tail', 'right');
    disp(['J mean rate test, 1 tail: ' num2str(x)]);
    x = signrank(mean_fits(:,2), 0, 'tail', 'right');
    disp(['P mean rate test, 1 tail: ' num2str(x)]);
    x = signrank(mean_fits(:,3), 0, 'tail', 'right');
    disp(['V mean rate test, 1 tail: ' num2str(x)]);
    
end