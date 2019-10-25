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
		  'ReturnVars',{''}, 'ArgsOnly',0, 'Cmds','', 'Map',0,'Smooth',1,'SameScale',0, ...
          'SIC',0,'Radii',0,'Details',0,'MinDur',0,'Filtered',0,'SortByRatio',0);
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
	n = Args.NumericArguments{1};
    if(Args.Map)
        disp('placeholder');
    end

	sdstr = get(obj,'SessionDirs');
	title(getDataOrder('ShortName','DirString',sdstr{n}))
else
	% plot all data
	n = get(obj,'Number');
    clf;
    if(Args.Map)
        if Args.Smooth
            maps = obj.data.maps_adsmooth;
            plot_all_maps(maps,Args);
        else
            maps = obj.data.maps_raw;
            plot_all_maps(maps,Args);
        end
    elseif(Args.SIC)
        histogram(obj.data.SICsh);
        max_count = max(histcounts(obj.data.SICsh));
        hold on;
        line([obj.data.SIC obj.data.SIC], [0 max_count]);
        hold off;
    elseif(Args.Details)
        detailed_plot(obj.data.detailed_fr,Args);
    elseif(Args.Radii)
        maps = obj.data.smoothed_size;
        plot_all_maps(maps,Args);
    end
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









function plot_all_maps(map,Args)

        % subplot('Position',[left bottom width height]);
        if Args.SameScale
            high = 0;
            for jj = 1:size(map,1)
                if max(max(map{jj})) > high
                    high = max(max(map{jj}));
                end           
            end
        end

        subplot('Position',[0.05 0.85 0.1 0.1]);
        if Args.SameScale
            imagesc(map{1}, [0 high]);
        else
            imagesc(map{1});
            colorbar();
        end
        title('1st cell');

        subplot('Position',[0.05 0.7 0.1 0.1]);
        if Args.SameScale
            imagesc(map{2}, [0 high]);
        else
            imagesc(map{2});
            colorbar();
        end
        title('2nd cell');

        subplot('Position',[0.20 0.65 0.3 0.3]);
        if Args.SameScale
            imagesc(map{3}, [0 high]);
        else
            imagesc(map{3});
            colorbar();
        end
        title('3rd cell');

        subplot('Position',[0.55 0.65 0.3 0.3]);
        if Args.SameScale
            imagesc(map{4}, [0 high]);
            colorbar();
        else
            imagesc(map{4});
            colorbar();
        end
        title('4th cell');

        subplot('Position',[0.05 0.5 0.9 0.09]);
        top_half = map{5};
        internal_max = max(max(top_half));
        cut = floor(size(top_half,2)/2);
        top_half = top_half(:,1:cut);
        bot_half = map{5};
        bot_half = bot_half(:,cut+1:end);
        if Args.SameScale
            imagesc(top_half, [0 high]);
        else
            imagesc(top_half, [0 internal_max]);
            colorbar();
        end
        title('5th cell (first then second half)');

        subplot('Position',[0.05 0.4 0.9 0.09]);
        if Args.SameScale
            imagesc(bot_half, [0 high]);
        else
            imagesc(bot_half, [0 internal_max]);
            colorbar();
        end

        subplot('Position',[0.025 0.25 0.45 0.1]);
        if Args.SameScale
            imagesc(map{6}, [0 high]);
        else
            imagesc(map{6});
            colorbar();
        end
        title('6th cell');

        subplot('Position',[0.525 0.25 0.45 0.1]);
        if Args.SameScale
            imagesc(map{7}, [0 high]);
        else
            imagesc(map{7});
            colorbar();
        end
        title('7th cell');

        subplot('Position',[0.025 0.1 0.45 0.1]);
        if Args.SameScale
            imagesc(map{8}, [0 high]);
        else
            imagesc(map{8});
            colorbar();
        end
        title('8th cell');

        subplot('Position',[0.525 0.1 0.45 0.1]);
        if Args.SameScale
            imagesc(map{9}, [0 high]);
        else
            imagesc(map{9});
            colorbar();
        end
        title('9th cell');
        
        
        
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
        disp(val);
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
