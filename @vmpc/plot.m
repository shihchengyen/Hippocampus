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
          'Map',0,'Smooth',1,'SIC',0,'Radii',0,'Details',0,'MinDur',0,'Filtered',0,'SortByRatio',0);
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
    if(Args.Errorbar)
    	if(obj.data.Args.UseMedian)
    		errorbar(1:size(obj.data.meanFRs,1),obj.data.meanFRs(:,n), ...
    			obj.data.semFRs(:,(n*2)-1),obj.data.semFRs(:,n*2),'x')
    	else
	        errorbar(obj.data.meanFRs(:,n),obj.data.semFRs(:,n),'.')
	    end
    elseif(Args.SIC)
        % get shuffled SIC
        shSIC = obj.data.SICsh(2:end,n);
        % get 95th percentile
        shSIC95 = prctile(shSIC,95);
        histogram(shSIC)
        hold on
        line(repmat(shSIC95,1,2),ylim,'Color','r')
        line(repmat(obj.data.SICsh(1,n),1,2),ylim,'Color','g')
        % line(repmat(obj.data.SIC(n),1,2),ylim,'Color','m')
        hold off
    elseif(Args.Shuffle)
        r = obj.data.SICsh(2:end,n);
        % get number of shuffles
        nShuffles = obj.data.Args.NumShuffles;
        ni = nShuffles/Args.ShuffleSteps;
        ra = zeros(ni,1);
        ii = 1:ni;
        for i = ii
            iend = i*Args.ShuffleSteps;
            ri = 1:iend;
            ra(i) = prctile(r(ri),95);
            if(iend(end)==1000)
                p1000 = ra(i);
            end
        end
        plot(ii*Args.ShuffleSteps,ra,'*-');
        hold on
        line(xlim,repmat(p1000,1,2),'Color','m')
        hold off
%         hold on
%         dra = 0.05 * mean(abs(diff(ra(1:4))));
%         line(repmat(xlim',1,2),repmat([ra(end)-dra ra(end)+dra],2,1))
%         hold off
        xlabel('Number of Shuffles')
        ylabel('95th percentile SIC')
    else  % if(Args.Errorbar)
        gSteps = obj.data.gridSteps;
        imagesc(reshape(obj.data.maps_raw(n,:),gSteps,gSteps));
        colorbar
    end  % if(Args.Errorbar)

	sdstr = get(obj,'SessionDirs');
	title(getDataOrder('ShortName','DirString',sdstr{n}))
else
	% plot all data
	n = get(obj,'Number');
    figure;
    if(Args.Map)
        if Args.Smooth
            maps = obj.data.maps_adsmooth;
            imagesc(reshape(maps,sqrt(length(maps)),sqrt(length(maps))));
            colorbar();
        else
            maps = obj.data.maps_raw;
            imagesc(reshape(maps,sqrt(length(maps)),sqrt(length(maps))));
            colorbar();
        end
    elseif(Args.SIC)
        histogram(obj.data.SICsh);
        max_count = max(histcounts(obj.data.SICsh));
        hold on;
        line([obj.data.SIC obj.data.SIC], [0 max_count]);
        hold off;
    elseif(Args.Details)
        detailed_plot(obj.data.detailed_fr,Args);
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
