function [obj, varargout] = plot(obj,varargin)
%@vmpc/plot Plot function for a vmpc object.
%   OBJ = plot(OBJ) creates an imagesc plot for the vmpc object. The
%   map shows the mean firing rates at each grid position.
%
%   InspectGUI(OBJ) creates an imagesc plot for a vmpc object
%       containing results from multiple cells. The GUI makes it easy to
%       plot the results for each cell.
%
%   InspectGUI(OBJ,'Errorbar') plots the results using an errorbar plot
%       instead of an imagesc plot.
%
%   InspectGUI(vpc,'addObjs',{vpc},'optArgs',{{},{'Errorbar'}},'SP',[2 1])
%       creates a figure with an imagesc plot on top of an errorbar plot.

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
		  'ReturnVars',{''}, 'ArgsOnly',0, 'Cmds','', 'Errorbar',0, ...
          'Shuffle',0, 'ShuffleSteps',100, 'NumSubPlots',4, ...
          'Map',0,'Smooth',1,'SIC',0,'Radii',0,'MinDur',0,'Filtered',1,...
          'SortByRatio',0,'Details',1,'RateBins',0,'MapOnly',0,'plotmap',0);
Args.flags = {'LabelsOff','ArgsOnly','Errorbar','SIC','Shuffle','MapOnly','plotmap'};
[Args,varargin2] = getOptArgs(varargin,Args);

if Args.plotmap
% Insert floor place map into larger 3D view setting
    
    mapL = obj.data.maps_adsm;
    ax = gca;
    mapLtemp = mapL;
    mapL = nan(1,5122);
    mapL(3:3+1600-1) = mapLtemp;
    mapG = flipud(reshape(mapLtemp, 40, 40)');
    mapGdummy = flipud(reshape(1:1600, 40, 40)');
       
    
    mapLdummy = 1:length(mapL);
    % Set up surf frame for plotting
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

    % ceiling follows floor mapping, top down view
    ceiling = flipud(reshape(mapL(1603:1603+1600-1), 40, 40)');
    ceilingdum = flipud(reshape(mapLdummy(1603:1603+1600-1), 40, 40)');

    % from top down, slit walls at bottom left corner, open outwards.
    % start from row closest to ground, rightwards, then climb rows
    walls = flipud(reshape(mapL(3203:3203+1280-1), 40*4, 8)');
    wallsdum = flipud(reshape(mapLdummy(3203:3203+1280-1), 40*4, 8)');

    % BL - bottom left, and so on, from top view, same slicing as walls
    % pillar width 8, height 5
    P1_BR = flipud(reshape(mapL(4483:4483+160-1), 8*4, 5)');
    P1_BRdum = flipud(reshape(mapLdummy(4483:4483+160-1), 8*4, 5)');
    P2_BL = flipud(reshape(mapL(4643:4643+160-1), 8*4, 5)');
    P2_BLdum = flipud(reshape(mapLdummy(4643:4643+160-1), 8*4, 5)');
    P3_TR = flipud(reshape(mapL(4803:4803+160-1), 8*4, 5)');
    P3_TRdum = flipud(reshape(mapLdummy(4803:4803+160-1), 8*4, 5)');
    P4_TL = flipud(reshape(mapL(4963:4963+160-1), 8*4, 5)');
    P4_TLdum = flipud(reshape(mapLdummy(4963:4963+160-1), 8*4, 5)');
    
 % Pad with NaNs for surf plots
    P1_BR = [P1_BR; nan(1,size(P1_BR,2))];
    P1_BR = [P1_BR nan(size(P1_BR,1),1)];

    P2_BL = [P2_BL; nan(1,size(P2_BL,2))];
    P2_BL = [P2_BL nan(size(P2_BL,1),1)];        

    P3_TR = [P3_TR; nan(1,size(P3_TR,2))];
    P3_TR = [P3_TR nan(size(P3_TR,1),1)];                

    P4_TL = [P4_TL; nan(1,size(P4_TL,2))];
    P4_TL = [P4_TL nan(size(P4_TL,1),1)];

    % Plot floor
    surf(floor_x, floor_y, floor_z, floor);
    alpha 1; shading flat;
    hold on;

    % Plot ceiling and walls
    surf(ceiling_x, ceiling_y, ceiling_z, ceiling);
    alpha 1; shading flat;
    surf(walls_x, walls_y, walls_z, walls);      
    alpha 1; shading flat;

    % Plot pillars
    surf(P1_x, P1_y, PX_z, P1_BR);
    alpha 1; shading flat;
    surf(P2_x, P2_y, PX_z, P2_BL);
    alpha 1; shading flat;
    surf(P3_x, P3_y, PX_z, P3_TR);
    alpha 1; shading flat;
    surf(P4_x, P4_y, PX_z, P4_TL);
    alpha 1; shading flat;

    axlim = max(abs([get(ax, 'xlim'), get(ax, 'ylim'),get(ax, 'zlim')])); % Make axes square around biggest value
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

% if user select 'ArgsOnly', return only Args structure for an empty object
if Args.ArgsOnly
    Args = rmfield (Args, 'ArgsOnly');
    varargout{1} = {'Args',Args};
    return;
end

if(~isempty(Args.NumericArguments))
	% plot one data set at a time
% 	n = get(gcf,'UserData');

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
        
    elseif(Args.Details)
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
        im0 = imagesc(reshape(map_choice,sqrt(length(map_choice)),sqrt(length(map_choice))), 'Tag','toppic');
        colorbar();  
        if ~Args.MapOnly     
            h1 = axes('Position',[0.1 0.1 0.8 0.3]);
            % what is the point of this plot?
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

% 	sdstr = get(obj,'SessionDirs');
% 	title(getDataOrder('ShortName','DirString',sdstr{n}))
else
	% plot all data
% 	n = get(obj,'Number');
    disp('placeholder');
end

if(~isempty(Args.Cmds))
    % save the current figure in case Args.Cmds switches to another figure
    h = gcf;
    eval(Args.Cmds)
    % switch back to previous figure
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
%         h0 = axes('Position',[0.3 0.5 0.4 0.4]);
        axes(h0);
        cla;
        set(h0,'Tag','top');
        if Args.Smooth
            map_choice = obj.data.maps_adsmooth(n,:);
        else
            map_choice = obj.data.maps_raw(n,:);
        end
        im0 = imagesc(reshape(map_choice,sqrt(length(map_choice)),sqrt(length(map_choice))), 'Tag','toppic');
        colorbar();
        
%         h1 = axes('Position',[0.1 0.1 0.8 0.3]);
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
        

function detailed_plot(details1,Args)

    unique_bins = unique(details1(1,:));
    if Args.Filtered
        checking_for_activity = details1(:,find(details1(2,:)>0));
        unique_bins = unique(checking_for_activity(1,:));
    end
    
    bin_limits = [0:5:25 max(details1(4,:))];
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
    

    for column = 1:min(10,size(binned_data,2))
        subplot(1,11,column);
        imagesc(binned_data(2:end,column));
        set(gca,'Ytick',1:size(binned_data(2:end,1),1)-1,'YTickLabel',[5:5:25]);
        colorbar();
        title(['grid ' num2str(unique_bins(column))]);
        xlabel({['zeros: ' num2str(binned_data(1,column))],['ratio: ' num2str(sum(binned_data(2:end,column))/binned_data(1,column))]});
    end
    
    if length(unique_bins) > 10
        ss = 1/floor(length(unique_bins)/10);
        uicontrol('Style', 'slider', 'Min', 0, 'Max', floor(length(unique_bins)/10), ...
            'Value', 0, 'SliderStep', [ss ss], 'Units','normalized', 'Position', [0.75 0.1 0.1 0.8], ...
            'Callback', {@react_to_slider, binned_data, unique_bins});
    end
    
    function react_to_slider(source, ~, binned_data, unique_bins)
        val = round(get(source, 'Value'));
        set(source, 'Value', val);
        counter = 1;
        for column = (10*val)+1:10*(val+1)
            sp_col = 11;
            subplot(1,sp_col,counter);
            if column > size(binned_data,2)
                data_to_clear = zeros(size(binned_data,1)-1,1);
                imagesc(data_to_clear);
                colorbar();
                title(' - ');
            else
                imagesc(binned_data(2:end,column));
                set(gca,'Ytick',1:size(binned_data(2:end,1),1)-1,'YTickLabel',[5:5:25]);
                colorbar();
                title(['grid ' num2str(unique_bins(column))]);
                xlabel({['zeros: ' num2str(binned_data(1,column))],['ratio: ' num2str(sum(binned_data(2:end,column))/binned_data(1,column))]});
            end
            counter = counter + 1;
        end        

        