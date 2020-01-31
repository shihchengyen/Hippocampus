function view_by_place_plot

    spiketrain = load('spiketrain.mat');
    spiketrain = spiketrain.timestamps ./ 1000; % in seconds
    
    pv = vmpv('auto');
    pv = pv.data.sessionTimeC;
    pv(:,4) = [0; diff(pv(:,1))];
    binned = histcounts(spiketrain, pv(:,1))';
    pv(:,5) = [binned; 0];
    
    pv(pv(:,2)==-1,:) = [];
    pv(pv(:,2)==0,:) = [];
    
    full_durations = NaN(5122,1600);
    full_spikes = NaN(5122,1600);
    
    for i = 1:1600
        
        subsample = [pv(pv(:,2)==i, [3 4 5])];
        
        % filling spikes
        subsample(subsample(:,3)==0,3) = nan;
        subsample(:,4) = circshift(subsample(:,2)~=0 ,-1);
        subsample(isnan(subsample(:,3)) & subsample(:,4), 3) = 0;
        subsample(:,4) = [];
        subsample(:,3) = fillmissing(subsample(:,3), 'next');
        
        % filling time
        subsample(subsample(:,2)==0,2) = nan;
        subsample(:,2) = fillmissing(subsample(:,2), 'previous');
        
        % padding with 5122 bin
        subsample = [subsample; [5122 0 0]];
       
        % remove bad view spots
        subsample(isnan(subsample(:,1)),:) = [];
        
        % sum durations
        full_durations(:,i) = accumarray(subsample(:,1), subsample(:,2));
        
%         % temp
%         full_durations(:,i) = accumarray(subsample(:,1), ones(size(subsample(:,1))));
        
%         % boundary checking
%         full_durations(end,i) = 1;
%         full_durations(5115,i) = 1;
%         full_durations(4482,i) = 1;
        
        % sum spikes
        full_spikes(:,i) = accumarray(subsample(:,1), subsample(:,3));
        
    end
    
    full_rate = full_spikes ./ full_durations;
    
    floor_x = repmat(0:40, 41, 1);
    floor_y = flipud(repmat([0:40]', 1, 41));
    floor_z = zeros(41,41);
    
    ceiling_x = floor_x;
    ceiling_y = floor_y;
    ceiling_z = 10.*ones(41,41);
    
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
    
    i = 1;    
    delta = 1;
    while 1==1
        
        placer = nan(1,1600);
        placer(i) = 1;
        template = nan(41,41);
        temp_c = template;
        temp_c(1:end-1,1:end-1) = flipud(reshape(placer(1:1600), 40, 40)');        
        place_c(:,:,1) = temp_c;        
        temp_c = template;
        temp_c(1:end-1,1:end-1) = 0.*flipud(reshape(placer(1:1600), 40, 40)');
        place_c(:,:,2) = temp_c;   
        temp_c = template;
        temp_c(1:end-1,1:end-1) = 0.*flipud(reshape(placer(1:1600), 40, 40)');     
        place_c(:,:,3) = temp_c;
        
        sampler = full_rate(:,i)';
%         sampler(sampler==0) = nan;
        
        if nansum(sampler) == 0 && 1==0
            continue;
        else
            disp(i);
        end
        
        % to view how mapping is done, replace sampler with 1:5122
        % floor starts from bottom left, going rightwards, climbing rows
        floor = flipud(reshape(sampler(3:3+1600-1), 40, 40)');
        % ceiling follows floor mapping, top down view
        ceiling = flipud(reshape(sampler(1603:1603+1600-1), 40, 40)');
        % from top down, slit walls at bottom left corner, open outwards.
        % start from row closest to ground, rightwards, then climb rows
        walls = flipud(reshape(sampler(3203:3203+1280-1), 40*4, 8)');
        % BL - bottom left, and so on, from top view, same slicing as walls
        % pillar width 8, height 5
        P1_BR = flipud(reshape(sampler(4483:4483+160-1), 8*4, 5)');
        P1_BR = [P1_BR; nan(1,size(P1_BR,2))];
        P1_BR = [P1_BR nan(size(P1_BR,1),1)];
        P2_BL = flipud(reshape(sampler(4643:4643+160-1), 8*4, 5)');
        P2_BL = [P2_BL; nan(1,size(P2_BL,2))];
        P2_BL = [P2_BL nan(size(P2_BL,1),1)];        
        P3_TR = flipud(reshape(sampler(4803:4803+160-1), 8*4, 5)');
        P3_TR = [P3_TR; nan(1,size(P3_TR,2))];
        P3_TR = [P3_TR nan(size(P3_TR,1),1)];                
        P4_TL = flipud(reshape(sampler(4963:4963+160-1), 8*4, 5)');
        P4_TL = [P4_TL; nan(1,size(P4_TL,2))];
        P4_TL = [P4_TL nan(size(P4_TL,1),1)];
        
        if i == 184
            save('debug.mat','P4_TL');
        end
        
        [caz, cel] = view;
        
        % floor
        surf(floor_x, floor_y, floor_z, floor);
        alpha 1; shading flat;
        hold on;
        
        % ceiling and walls
        surf(ceiling_x, ceiling_y, ceiling_z, ceiling);
        alpha 1; shading flat;
        surf(walls_x, walls_y, walls_z, walls);      
        alpha 1; shading flat;
        
        disp(sum(sum(find(ceiling==0))) + sum(sum(find(floor==0))) + sum(sum(find(P4_TL==0))));
        
        % pillars
        surf(P1_x, P1_y, PX_z, P1_BR);
        alpha 1; shading flat;
        surf(P2_x, P2_y, PX_z, P2_BL);
        alpha 1; shading flat;
        surf(P3_x, P3_y, PX_z, P3_TR);
        alpha 1; shading flat;
        surf(P4_x, P4_y, PX_z, P4_TL);
        alpha 1; shading flat;    
        
        % marking out place
        place_plot = surf(floor_x, floor_y, floor_z, nan(41,41));
        place_plot.CData = place_c;
        shading flat;
        
        % marking outlines
        plot3([8 16 nan 24 32], [8 8 nan 8 8], [5 5 nan 5 5], 'k');
        plot3([8 16 nan 24 32], [8 8 nan 8 8], [0 0 nan 0 0], 'k');
        plot3([8 16 nan 24 32], [16 16 nan 16 16], [5 5 nan 5 5], 'k');
        plot3([8 16 nan 24 32], [16 16 nan 16 16], [0 0 nan 0 0], 'k');
        plot3([8 16 nan 24 32], [24 24 nan 24 24], [5 5 nan 5 5], 'k');
        plot3([8 16 nan 24 32], [24 24 nan 24 24], [0 0 nan 0 0], 'k');
        plot3([8 16 nan 24 32], [32 32 nan 32 32], [5 5 nan 5 5], 'k');
        plot3([8 16 nan 24 32], [32 32 nan 32 32], [0 0 nan 0 0], 'k');
        plot3([8 8 nan 8 8], [8 16 nan 24 32], [5 5 nan 5 5], 'k');
        plot3([8 8 nan 8 8], [8 16 nan 24 32], [0 0 nan 0 0], 'k');
        plot3([16 16 nan 16 16], [8 16 nan 24 32], [5 5 nan 5 5], 'k');
        plot3([16 16 nan 16 16], [8 16 nan 24 32], [0 0 nan 0 0], 'k');
        plot3([24 24 nan 24 24], [8 16 nan 24 32], [5 5 nan 5 5], 'k');
        plot3([24 24 nan 24 24], [8 16 nan 24 32], [0 0 nan 0 0], 'k');
        plot3([32 32 nan 32 32], [8 16 nan 24 32], [5 5 nan 5 5], 'k');
        plot3([32 32 nan 32 32], [8 16 nan 24 32], [0 0 nan 0 0], 'k');  
        plot3([8 8 nan 8 8 nan 8 8 nan 8 8], [8 8 nan 16 16 nan 24 24 nan 32 32], [0 5 nan 0 5 nan 0 5 nan 0 5], 'k');
        plot3([16 16 nan 16 16 nan 16 16 nan 16 16], [8 8 nan 16 16 nan 24 24 nan 32 32], [0 5 nan 0 5 nan 0 5 nan 0 5], 'k');
        plot3([24 24 nan 24 24 nan 24 24 nan 24 24], [8 8 nan 16 16 nan 24 24 nan 32 32], [0 5 nan 0 5 nan 0 5 nan 0 5], 'k');
        plot3([32 32 nan 32 32 nan 32 32 nan 32 32], [8 8 nan 16 16 nan 24 24 nan 32 32], [0 5 nan 0 5 nan 0 5 nan 0 5], 'k');
        
        hold off; colormap(gca, copper); colorbar;
        view(gca, [caz cel]);
        
        input_num = input('continue, 0 to reverse, 2 4 6 8 to traverse\n');
        if isempty(input_num)
            i = i + delta;
        elseif input_num == 0
            delta = -1*delta;
            i = i + delta;
        elseif input_num == 52
            i = i - 5*40;            
        elseif input_num == 2
            i = i - 40;
        elseif input_num == 8
            i = i + 40;
        elseif input_num == 58
            i = i + 5*40;
        elseif input_num == 4
            i = i - 1;
        elseif input_num == 54
            i = i - 5;            
        elseif input_num == 6
            i = i + 1;
        elseif input_num == 56
            i = i + 5;            
        end
        if i > 1600
            i = 1600;
        end
        if i < 1
            i = 1;
        end
        
    end
        
        
end



