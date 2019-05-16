function [obj, varargout] = plot(obj,varargin)
%@vmplacecell/plot Plot function for a vmplacecell object.
%   OBJ = plot(OBJ) creates an imagesc plot for the vmplacecell object. The
%   map shows the mean firing rates at each grid position.
%
%   InspectGUI(OBJ) creates an imagesc plot for a vmplacecell object
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
          'SIC',0,'Shuffle',0, 'ShuffleSteps',100, 'NumSubPlots',4);
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
        imagesc(reshape(obj.data.meanFRs(:,n),gSteps,gSteps));
        colorbar
    end  % if(Args.Errorbar)

	sdstr = get(obj,'SessionDirs');
	title(getDataOrder('ShortName','DirString',sdstr{n}))
else
	% plot all data
	n = get(obj,'Number');
    if(Args.SIC)
        % get data SIC
        dSIC = obj.data.SICsh(1,:);
        sSIC = obj.data.SICsh(2:end,:);
        sSIC95 = prctile(sSIC,95);
        % plot SIC versus shuffle
        nsp = Args.NumSubPlots;
        % get number of cells per subplot
        npsubplot = ceil(n/nsp);
        % generate indices
        npi = 1:npsubplot;
        for i = 1:nsp
            subplot(nsp,1,i)
            dindices = npi+((i-1)*npsubplot);
            plot(dindices,dSIC(dindices),'b*-')
            hold on
            plot(dindices,sSIC95(dindices),'ro-')
            hold off
            xlim([dindices(1)-1 dindices(end)+1])
        end
        xlabel('Cell Number')
        ylabel('SIC')
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
