function [obj, varargout] = plot(obj,varargin)
%@vmpv/plot Plot function for a vmpv object.
%   OBJ = plot(OBJ) creates an imagesc plot for the vmpv object. The
%   map shows the mean firing rates at each grid position.
%
%   InspectGUI(OBJ) creates an imagesc plot for a vmpv object
%       containing results from multiple cells. The GUI makes it easy to
%       plot the results for each cell.
%
%   InspectGUI(OBJ,'Errorbar') plots the results using an errorbar plot
%       instead of an imagesc plot.
%
%   InspectGUI(vpc,'addObjs',{vpc},'optArgs',{{},{'Errorbar'}},'SP',[2 1])
%       creates a figure with an imagesc plot on top of an errorbar plot.

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
		  'ReturnVars',{''}, 'ArgsOnly',0, 'Cmds','');
Args.flags = {'LabelsOff','ArgsOnly','Errorbar','SIC','Shuffle'};
[Args,varargin2] = getOptArgs(varargin,Args);

% if user select 'ArgsOnly', return only Args structure for an empty object
if Args.ArgsOnly
    Args = rmfield (Args, 'ArgsOnly');
    varargout{1} = {'Args',Args};
    return;
end

if(~isempty(Args.NumericArguments))
	% plot one data set at a time
% 	n = get(gcf,'UserData');

    po = findobj(gcf,'String','Plot Options');
    set(po, 'Visible', 'off');

    % plotting not yet implemented

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

        