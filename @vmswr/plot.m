function [obj, varargout] = plot(obj,varargin)
%@vmswr/plot Plot function for vmswr object.
%   OBJ = plot(OBJ) creates a raster plot of the neuronal
%   response.

Args = struct('LabelsOff',0,'TitleOff', 0,'Color','b', ...
            'GroupPlots',1,'GroupPlotIndex',1,...
            'PreTrial',500,'Rms',0,'Cmds','','DisplayTrials',0,...
            'PlotAllData',0,'RemoveLineNoise',[],...
            'ReturnVars',{''}, 'ArgsOnly',0);
Args.flags = {'LabelsOff','TitleOff','ArgsOnly', ...
            'PlotAllData','Rms','DisplayTrials'};
[Args,varargin2] = getOptArgs(varargin,Args);

% if user select 'ArgsOnly', return only Args structure for an empty object
if Args.ArgsOnly
    Args = rmfield (Args, 'ArgsOnly');
    varargout{1} = {'Args',Args};
    return;
end

sRate = obj.data.analogInfo.SampleRate;
swrInfo = obj.data.analogRmsInfo.Swr;

% for title
sdstr = get(obj,'SessionDirs');
[a,b] = fileparts(sdstr{1});
title(b)
sp = {' '}; %space

if(~isempty(Args.NumericArguments))
	% plot one data set at a time
	n = Args.NumericArguments{1};

    tIdx = obj.data.trialIndices(n,:);
    swrMean = obj.data.analogRmsInfo.Mean;
    swrStd = obj.data.analogRmsInfo.Std;
    
    if (~Args.PlotAllData)
        if(length(tIdx)==2)
            OldMarkerFormat = 1;
        else
            OldMarkerFormat = 0;
        end

        if(OldMarkerFormat)
            idx = (tIdx(1)-(Args.PreTrial/1000*sRate)):tIdx(2);
        else
            idx = (tIdx(1)-(Args.PreTrial/1000*sRate)):tIdx(3);
        end
               
    else
        idx = 1:length(obj.data.analogData);
    end
    
    tStart = idx(1); % for 
    tEnd = idx(end);
    swridx = unique(roundn(idx./5,0));

    if(Args.Rms)
        sc = 1000/200; %scaling for graphs
        idx=swridx;
        if (Args.PlotAllData)
            data = obj.data.analogRmsData';
        else
            data = obj.data.analogRmsData';
            data = data(idx,:);
        end
        numPlots = floor(size(data,1)/200);
    else
        sc = 1000/1000;
        if (Args.PlotAllData)
            data = obj.data.analogData;
        else
            data = obj.data.analogData(idx,:);
        end
        
        if(~isempty(Args.RemoveLineNoise))
			data = nptRemoveLineNoise(data,Args.RemoveLineNoise,sRate);
		end
    end

    if(~isempty(Args.RemoveLineNoise))
            data = nptRemoveLineNoise(data,Args.RemoveLineNoise,sRate);
    end

    % set axis
    y0 = min(data);
    y1 = max(data);

    di = idx(1);
    df = idx(end);

    % if (~isempty(gco))
    %     delete(gca)
    % end

    hold on

    for i = 1:length(swrInfo)
        if swrInfo(i,3) >= swridx(1) && swrInfo(i,2) <= swridx(end) % detects if swr event is in the graph
    %             si = 0; %swr event
            if swrInfo(i,3) >= swridx(1) && swrInfo(i,2) < swridx(1)
                x0 = swridx(1);
                x1 = swrInfo(i,3)+0.5;
    %                 if swrInfo(i,1) >= swridx(1)
    %                     si = swrInfo(i,1);
    %                 end
            elseif swrInfo(i,2) <= swridx(end) && swrInfo(i,3) > swridx(end)
                x0 = swrInfo(i,2)-0.5;
                x1 = swridx(end);
    %                 if swrInfo(i,1) <= swridx(end)
    %                     si = swrInfo(i,1);
    %                 end
            else
                x0 = swrInfo(i,2)-0.5;
                x1 = swrInfo(i,3)+0.5;
    %                 si = swrInfo(i,1);
            end
    %             si = (si-swridx(1)).*sc;
            patch(([x0 x1 x1 x0]).*5, [y0*[1 1] y1*[1 1]], [1 0.1 0.1], 'LineStyle','none')
    %             if (si ~= 0)
    %                 plot(si, data(si/sc),'rx')
    %             end
        end
    end

    plot((di:df)*sc, data,'Color',Args.Color);

    ylim([y0 y1])
    xlim([di*sc df*sc])

    % plotting std away from mean
    plot(xlim,[swrMean+3*swrStd swrMean+3*swrStd],'g');
    plot(xlim,[swrMean+swrStd swrMean+swrStd],'g'); 
    text(di*sc, swrMean+3*swrStd, 'threshold = 3\sigma');
    text(di*sc, swrMean+swrStd, '1\sigma');
    % p = gco;
    hold off
    
    % texts for title, yaxis and xaxis 
    tt = strcat('LFP of',sp,b,sp,'from',sp,num2str(tStart),'ms to',sp,num2str(tEnd),'ms');
    ty = 'LFP';
    tx = 'time (ms)';
else
  
    tIdx = obj.data.trialIndices;
    s = size(tIdx);
    tIdx(:,1) = (tIdx(:,1)-(Args.PreTrial/1000*sRate))./5;
    
    if(s(2)==2)
        OldMarkerFormat = 1;
    else
        OldMarkerFormat = 0;
    end
    
    if(OldMarkerFormat)
        tIdx(:,2) = (tIdx(:,2))./5;
    else
        tIdx(:,2) = (tIdx(:,3))./5;
    end
    
    % count SPW Event in each trial 
    for i = 1:length(tIdx)
        tIdx(i,3) = sum(swrInfo(:,1)>tIdx(i,1) & swrInfo(:,1)<tIdx(i,2));
    end
    
    histogram(tIdx(:,3))
    
    % texts for title, yaxis and xaxis 
    tt = strcat('Histogram of SPW Events in',sp,b,sp,'in each trial');
    ty = 'Trial';
    tx = 'SPW Event (count)';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @vmswr/PLOT takes 'LabelsOff' and 'TitleOff'
Args.LabelsOff = 0;

if(~Args.LabelsOff)
    xlabel(tx)
    ylabel(ty)
end

if(~Args.TitleOff)
    title(tt)
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(~isempty(Args.Cmds))
    % save the current figure in case Args.Cmds switches to another figure
%     h = gcf;
    % save current directory
%     cwd = pwd;
%     cd(sdstr{1})
    eval(Args.Cmds)
    % switch back to previous figure
%     figure(h);
%     cd(cwd)
end

if(Args.DisplayTrials)
    vr = vmswr('auto');
    InspectGUI(vr,'addObjs',{vr},'optArgs',{{},{'Rms'}},'SP',[2 1])
end

% hidata plot overwrites the x and y label
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