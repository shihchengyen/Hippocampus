function [obj, varargout] = plot(obj,varargin)
%@vmpc/plot Plot function for a vmpc object.
%   Use InspectGUI(OBJ) to create an interactive plot interface

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
		  'ReturnVars',{''}, 'ArgsOnly',0, 'Cmds','', 'Errorbar',0, ...
          'Shuffle',0, 'ShuffleSteps',100, 'NumSubPlots',4, ...
          'Map',0,'Smooth',1,'SIC',0,'Radii',0,'MinDur',0,'Filtered',1,...
          'SortByRatio',0,'Details',1,'RateBins',0,'MapOnly',0,'plotmap',0, ...
          'Occupancy', 0, 'outline',0, 'outline_threshold',0.75, 'peak_detection_pct',0.85);
Args.flags = {'LabelsOff','ArgsOnly','Errorbar','SIC','Shuffle','MapOnly','plotmap', 'Occupancy', 'outline'};
[Args,varargin2] = getOptArgs(varargin,Args);

figure; 

if Args.plotmap
    mapL = obj.data.maps_adsm;
    ax = gca;
    mapLtemp = mapL;
    mapL = nan(1,5122);
    mapL(3:3+1600-1) = mapLtemp;
    mapG = flipud(reshape(mapLtemp, 40, 40)');
    mapGdummy = flipud(reshape(1:1600, 40, 40)');
    mapLdummy = 1:length(mapL);
    floor_x = repmat(0:40, 41, 1);
    floor_y = flipud(repmat([0:40]', 1, 41));
    floor_z = zeros(41,41);
    ceiling_x = floor_x;
    ceiling_y = floor_y;
    ceiling_z = 40.*ones(41,41);
    walls_x = repmat([0.*ones(1,40) 0:39 40.*ones(1,40) 40:-1:0], 9, 1);
    walls_y = repmat([0:39 40.*ones(1,40) 40:-1:1 0.*ones(1,41)], 9, 1);
    walls_z = repmat([24:-1:16]', 1, 40*4 + 1);
    P1_x = repmat([24.*ones(1,8) 24:31 32.*ones(1,8) 32:-1:24], 6, 1);
    P1_y = repmat([8:15 16.*ones(1,8) 16:-1:9 8.*ones(1,9)], 6, 1);
    PX_z = repmat([21:-1:16]', 1, 8*4 + 1);
    P2_x = repmat([8.*ones(1,8) 8:15 16.*ones(1,8) 16:-1:8], 6, 1);
    P2_y = P1_y;
    P3_x = P1_x;
    P3_y = repmat([24:31 32.*ones(1,8) 32:-1:25 24.*ones(1,9)], 6, 1);
    P4_x = P2_x;
    P4_y = P3_y;
    floor = flipud(reshape(mapL(3:3+1600-1), 40, 40)');
    floordum = flipud(reshape(mapLdummy(3:3+1600-1), 40, 40)');
    ceiling = flipud(reshape(mapL(1603:1603+1600-1), 40, 40)');
    ceilingdum = flipud(reshape(mapLdummy(1603:1603+1600-1), 40, 40)');
    walls = flipud(reshape(mapL(3203:3203+1280-1), 40*4, 8)');
    wallsdum = flipud(reshape(mapLdummy(3203:3203+1280-1), 40*4, 8)');
    P1_BR = flipud(reshape(mapL(4483:4483+160-1), 8*4, 5)');
    P1_BRdum = flipud(reshape(mapLdummy(4483:4483+160-1), 8*4, 5)');
    P2_BL = flipud(reshape(mapL(4643:4643+160-1), 8*4, 5)');
    P2_BLdum = flipud(reshape(mapLdummy(4643:4643+160-1), 8*4, 5)');
    P3_TR = flipud(reshape(mapL(4803:4803+160-1), 8*4, 5)');
    P3_TRdum = flipud(reshape(mapLdummy(4803:4803+160-1), 8*4, 5)');
    P4_TL = flipud(reshape(mapL(4963:4963+160-1), 8*4, 5)');
    P4_TLdum = flipud(reshape(mapLdummy(4963:4963+160-1), 8*4, 5)');
    P1_BR = [P1_BR; nan(1,size(P1_BR,2))];
    P1_BR = [P1_BR nan(size(P1_BR,1),1)];
    P2_BL = [P2_BL; nan(1,size(P2_BL,2))];
    P2_BL = [P2_BL nan(size(P2_BL,1),1)];        
    P3_TR = [P3_TR; nan(1,size(P3_TR,2))];
    P3_TR = [P3_TR nan(size(P3_TR,1),1)];                
    P4_TL = [P4_TL; nan(1,size(P4_TL,2))];
    P4_TL = [P4_TL nan(size(P4_TL,1),1)];
    surf(floor_x, floor_y, floor_z, floor);
    alpha 1; shading flat;
    hold on;
    surf(ceiling_x, ceiling_y, ceiling_z, ceiling);
    alpha 1; shading flat;
    surf(walls_x, walls_y, walls_z, walls);      
    alpha 1; shading flat;
    surf(P1_x, P1_y, PX_z, P1_BR);
    alpha 1; shading flat;
    surf(P2_x, P2_y, PX_z, P2_BL);
    alpha 1; shading flat;
    surf(P3_x, P3_y, PX_z, P3_TR);
    alpha 1; shading flat;
    surf(P4_x, P4_y, PX_z, P4_TL);
    alpha 1; shading flat;
    axlim = max(abs([get(ax, 'xlim'), get(ax, 'ylim'),get(ax, 'zlim')]));
    axis(ax,[0 axlim 0 axlim 0 axlim]);  
    colormap jet;
    alpha 1; shading flat; 
    view(-35,20);
    axis(ax,'tight');
    w = ax.OuterPosition(3)*0.02;
    x = ax.OuterPosition(1) + ax.OuterPosition(3);
    y = ax.OuterPosition(2) + ax.OuterPosition(4)/4;
    h = ax.OuterPosition(4)/2;
    c = colorbar(ax,'Position',[x y w h]);
end

if Args.ArgsOnly
    Args = rmfield (Args, 'ArgsOnly');
    varargout{1} = {'Args',Args};
    return;
end

if Args.Occupancy
    if Args.Smooth
        dur_map = obj.data.dur_adsm;
    else
        dur_map = obj.data.dur_raw;
    end
    if length(dur_map) > 1600
        dur_floor = dur_map(3:1602);
    else
        dur_floor = dur_map;
    end
    imagesc(flipud(reshape(dur_floor, 40, 40)'));
    colorbar;
    title('Occupancy Map (Duration) - Floor Only');
    xlabel('X Position');
    ylabel('Y Position');
    axis square;
    colormap jet;
end

if(~isempty(Args.NumericArguments))
    colormap jet;
    po = findobj(gcf,'String','Plot Options');
    set(po, 'Visible', 'off');
    n = Args.NumericArguments{1};
    text_field = findobj(gcf,'Tag','StaticText1');
    set(text_field,'UserData',n);
    disp_number = findobj(gcf,'Tag','EditText1');
    set(disp_number,'String',num2str(n));

    if(Args.SIC)
        histogram(obj.data.SICsh(:,n));
        max_count = max(histcounts(obj.data.SICsh(:,n)));
        hold on;
        line([obj.data.SIC(n,1) obj.data.SIC(n,1)], [0 max_count], 'color','red');
        hold off;
        title(obj.data.origin{n});
        next_handle = findobj(gcf,'String','Next');
        prev_handle = findobj(gcf,'String','Previous');
        set(next_handle,'Callback',{@forwardcallback, length(obj.data.origin), Args, gcf, obj, 'SIC'});
        set(prev_handle,'Callback',{@backcallback, Args, gcf, obj, 'SIC'});         
        
    elseif(Args.Details && ~Args.Occupancy)
        set(gca,'visible','off');
        
        if Args.MapOnly
            h0 = axes('Position',[0.1 0.1 0.8 0.8]);
        else
            h0 = axes('Position',[0.3 0.5 0.4 0.4]);
        end
        set(h0,'Tag','top');
        
        if Args.Smooth
            map_choice = obj.data.maps_adsm(n,:);
        else
            map_choice = obj.data.maps_raw(n,:);
        end
        
        grid_size = sqrt(length(map_choice));
        map_2d = reshape(map_choice, grid_size, grid_size); 
        
        if Args.outline
            % Get outline parameters
            threshold_percentile = Args.outline_threshold;
            peak_pct = Args.peak_detection_pct;
            
            % Call get_outline with both parameters
            [region_cells, peak_indices, peak_values, ~, ~] = get_outline(obj, ...
                'outline_threshold', threshold_percentile, ...
                'peak_detection_pct', peak_pct, ...
                'Smooth', Args.Smooth, ...
                'CellIndex', n);
            
            if ~isempty(region_cells)
                % Create RGB image for visualization
                rgb_image = zeros(grid_size, grid_size, 3);
                
                % Initialize background
                for r = 1:grid_size
                    for c = 1:grid_size
                        if isnan(map_2d(r,c))
                            rgb_image(r,c,:) = [1 1 1];  % White for non-maze
                        else
                            rgb_image(r,c,:) = [0.7 0.7 0.7];  % Gray for maze
                        end
                    end
                end
                
                % Color all regions in RED
                region_color = [1.0, 0.0, 0.0];
                num_regions = length(region_cells);
                for reg = 1:num_regions
                    region_idx = region_cells{reg};
                    [rows, cols] = ind2sub([grid_size, grid_size], region_idx);
                    for i = 1:length(rows)
                        rgb_image(rows(i), cols(i), :) = region_color;
                    end
                end
                
                % Use imagesc instead of image for proper scaling
                im0 = imagesc(rgb_image);
                set(im0, 'Tag', 'toppic');
                axis square;
                
                hold on;
                % Mark all peaks with white stars
                for reg = 1:num_regions
                    [peak_row, peak_col] = ind2sub([grid_size, grid_size], peak_indices(reg));
                    plot(peak_col, peak_row, 'w*', 'MarkerSize', 10, 'LineWidth', 2);
                    text(peak_col, peak_row-1, sprintf('%d', reg), ...
                        'Color', 'white', 'FontSize', 10, 'FontWeight', 'bold', ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
                end
                hold off;
                
                % Title
                total_pixels = sum(cellfun(@length, region_cells));
                title(sprintf('Cell %d | %d Regions | Display: %.0f%% | Peak Det: %.0f%%', ...
                    n, num_regions, threshold_percentile*100, peak_pct*100));
                
            else
                % No regions found
                im0 = imagesc(map_2d, 'Tag', 'toppic');
                axis square;
                colormap jet;
                colorbar();
            end
        else
            % Normal plot without outline
            im0 = imagesc(map_2d, 'Tag', 'toppic');
            axis square;
            colorbar();
        end
        
        if ~Args.MapOnly && isfield(obj.data, 'detailed_fr')  
            h1 = axes('Position',[0.1 0.1 0.8 0.3]);
            details1 = obj.data.detailed_fr{n,1};
            unique_bins = unique(details1(1,:));
            if Args.Filtered
                checking_for_activity = details1(:,find(details1(2,:)>0));
                unique_bins = unique(checking_for_activity(1,:));
            end
            bin_limits = [0:5:25 max(details1(4,:))];
            if Args.RateBins ~= 0
                if length(Args.RateBins) == 1
                    if Args.RateBins > 0
                        bin_limits = [0:Args.RateBins:25 max(details1(4,:))];
                    else
                        bin_limits = [exp(0.1.*[0:25])-1 max(details1(4,:))];
                    end 
                end
            end
            binned_data = NaN(length(bin_limits),length(unique_bins));
            for col = 1:length(unique_bins)
                subset_arr = details1(3:4,find(details1(1,:)==unique_bins(col)));
                if Args.MinDur ~= 0
                    subset_arr = subset_arr(:,find(subset_arr(1,:)>Args.MinDur));
                end
                zero_count = sum(subset_arr(2,:)==0);
                binned_temp = histcounts(subset_arr(2,:), bin_limits);
                binned_temp(1) = binned_temp(1) - zero_count;
                binned_data(1:length(binned_temp)+1,col) = [zero_count binned_temp];
            end
            if Args.SortByRatio
                ratio = sum(binned_data(2:end,:),1)./binned_data(1,:);
                binned_data_temp = [ratio; binned_data; unique_bins];
                binned_data_temp = sortrows(binned_data_temp.',1).';
                binned_data = fliplr(binned_data_temp(2:end-1,:));
                unique_bins = fliplr(binned_data_temp(end,:));
            end
            
            im1 = imagesc(binned_data(2:end,:), 'Tag','botpic');
            set(gca,'YTick',1:length(bin_limits)-2,'YTickLabel',bin_limits(2:end-1));
            title(obj.data.origin{n});
            colorbar();
            
            next_handle = findobj(gcf,'String','Next');
            prev_handle = findobj(gcf,'String','Previous');
            set(next_handle,'Callback',{@forwardcallback, length(obj.data.origin), Args, gcf, obj, h0, h1});
            set(prev_handle,'Callback',{@backcallback, Args, gcf, obj, h0, h1});            
            set(gcf,'WindowButtonMotionFcn',{@hovercallback,unique_bins, bin_limits,sqrt(length(map_choice)),h0,h1,binned_data,map_choice,im0,im1});
        end
    end
else
    disp('placeholder');
end

if(~isempty(Args.Cmds))
    h = gcf;
    eval(Args.Cmds)
    figure(h);
end

RR = eval('Args.ReturnVars');
lRR = length(RR);
if(lRR>0)
    for i=1:lRR
        RR1{i}=eval(RR{i});
    end 
    varargout = getReturnVal(Args.ReturnVars, RR1);
else
    varargout = {};
end

% Callback functions
function hovercallback(source, ~, unique_bins,bin_limits,dim,h0,h1,binned_data,full_map,im0,im1)
    hAxes = hittest(gcf);
    hover_loc = get(hAxes, 'Tag');
    if strcmpi(hover_loc,'toppic')
        cpt = get(h0,'CurrentPoint');
        ad = ones(1,dim*dim);
        ad(1,(floor(cpt(1,1))-1)*dim + floor(cpt(1,2))) = 0;
        hAxes.AlphaData = reshape(ad, dim, dim);  
        ad = ones(size(im1.CData));
        index = floor(cpt(1,1))*dim + floor(cpt(1,2));
        loc = find(unique_bins==index);
        if ~isempty(loc)
            ad(:,loc) = 0;
            im1.AlphaData = ad;
        end
    elseif strcmpi(hover_loc,'botpic')
        cpt = get(h1,'CurrentPoint');
        ad = ones(size(hAxes.CData));
        ad(:,round(cpt(1,1))) = 0;
        hAxes.AlphaData = ad;
        ad = ones(1,dim*dim);
        ad(1,unique_bins(round(cpt(1,1)))) = 0;
        im0.AlphaData = reshape(ad, dim, dim);          
    end

function [Args, gcf, obj] = forwardcallback(source, ~, limit, Args, gcf, obj, h0, h1)
    text_field = findobj(gcf,'Tag','StaticText1');
    disp_number = findobj(gcf,'Tag','EditText1');
    n = get(text_field,'UserData');
    if n < limit
        n = n + 1;
        set(text_field, 'UserData', n);
        set(disp_number,'String',num2str(n));
    end
    if ~strcmpi(h0, 'SIC')
        replot(gcf, Args, obj, h0, h1);
    else
        histogram(obj.data.SICsh(:,n));
        max_count = max(histcounts(obj.data.SICsh(:,n)));
        hold on;
        line([obj.data.SIC(n,1) obj.data.SIC(n,1)], [0 max_count],'color','red');
        hold off;
        title(obj.data.origin{n});
    end

function backcallback(source, ~, Args, gcf, obj, h0, h1)
    text_field = findobj(gcf,'Tag','StaticText1');
    disp_number = findobj(gcf,'Tag','EditText1');
    n = get(text_field,'UserData');
    if n > 1
        n = n - 1;
        set(text_field, 'UserData', n);
        set(disp_number,'String',num2str(n));
    end
    if ~strcmpi(h0, 'SIC')
        replot(gcf, Args, obj, h0, h1);
    else
        histogram(obj.data.SICsh(:,n));
        max_count = max(histcounts(obj.data.SICsh(:,n)));
        hold on;
        line([obj.data.SIC(n,1) obj.data.SIC(n,1)], [0 max_count],'color','red');
        hold off;
        title(obj.data.origin{n});
    end

function replot(gcf, Args, obj, h0, h1)
    text_field = findobj(gcf,'Tag','StaticText1');
    n = get(text_field,'UserData');
    if(Args.SIC)
        histogram(obj.data.SICsh(:,n));
        max_count = max(histcounts(obj.data.SICsh(:,n)));
        hold on;
        line([obj.data.SIC(n,1) obj.data.SIC(n,1)], [0 max_count]);
        hold off;
    elseif(Args.Details)
        set(gca,'visible','off');
        axes(h0);
        cla;
        set(h0,'Tag','top');
        if Args.Smooth
            map_choice = obj.data.maps_adsm(n,:);
        else
            map_choice = obj.data.maps_raw(n,:);
        end
        im0 = imagesc(reshape(map_choice,sqrt(length(map_choice)),sqrt(length(map_choice))), 'Tag','toppic');
        colorbar();
        axes(h1);
        cla;
        details1 = obj.data.detailed_fr{n,1};
        unique_bins = unique(details1(1,:));
        if Args.Filtered
            checking_for_activity = details1(:,find(details1(2,:)>0));
            unique_bins = unique(checking_for_activity(1,:));
        end
        bin_limits = [0:5:25 max(details1(4,:))];
        if Args.RateBins ~= 0
            if length(Args.RateBins) == 1
                if Args.RateBins > 0
                    bin_limits = [0:Args.RateBins:25 max(details1(4,:))];
                else
                    bin_limits = [exp(0.1.*[0:25])-1 max(details1(4,:))];
                end
            end
        end
        binned_data = NaN(length(bin_limits),length(unique_bins));
        for col = 1:length(unique_bins)
            subset_arr = details1(3:4,find(details1(1,:)==unique_bins(col)));
            if Args.MinDur ~= 0
                subset_arr = subset_arr(:,find(subset_arr(1,:)>Args.MinDur));
            end
            zero_count = sum(subset_arr(2,:)==0);
            binned_temp = histcounts(subset_arr(2,:), bin_limits);
            binned_temp(1) = binned_temp(1) - zero_count;
            binned_data(1:length(binned_temp)+1,col) = [zero_count binned_temp];
        end
        if Args.SortByRatio
            ratio = sum(binned_data(2:end,:),1)./binned_data(1,:);
            binned_data_temp = [ratio; binned_data; unique_bins];
            binned_data_temp = sortrows(binned_data_temp.',1).';
            binned_data = fliplr(binned_data_temp(2:end-1,:));
            unique_bins = fliplr(binned_data_temp(end,:));
        end
        im1 = imagesc(binned_data(2:end,:), 'Tag','botpic');
        set(gca,'YTick',1:length(bin_limits)-2,'YTickLabel',bin_limits(2:end-1));
        title(obj.data.origin{n});
        colorbar();
        set(gcf,'WindowButtonMotionFcn',{@hovercallback,unique_bins, bin_limits,sqrt(length(map_choice)),h0,h1,binned_data,map_choice,im0,im1});
    end