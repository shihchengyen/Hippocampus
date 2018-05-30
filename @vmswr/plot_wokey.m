function [obj, varargout] = plot(obj,varargin)
%@vmswr/plot Plot function for vmswr object.
%   OBJ = plot(OBJ) creates a raster plot of the neuronal
%   response.

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
            'PreTrial',500,'Rms',0, ...
            'PlotAllData',0,'RemoveLineNoise',[],'TitleOff', 0, ...
            'ReturnVars',{''}, 'ArgsOnly',0);
Args.flags = {'LabelsOff','ArgsOnly', ...
            'PlotAllData'};
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
else
	% plot all data
	n = 1;
    fprintf('work in progress');
end

% find indices for n-th trial
tIdx = obj.data.trialIndices(n,:);
sRate = obj.data.analogInfo.SampleRate;
swrInfo = obj.data.analogRmsInfo.Swr;
swrMean = obj.data.analogRmsInfo.Mean;
swrStd = obj.data.analogRmsInfo.Std;

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

if(Args.Rms)
    idx=unique(roundn(idx./5,0));
    if (Args.PlotAllData)
        data = obj.data.analogRmsData';
    else
        data = obj.data.analogRmsData';
        data = data(idx,:);
    end
    sc = 200; %scaling for graphs
    numPlots = floor(size(data,1)/sc);
else
    if (Args.PlotAllData)
        data = obj.data.analogData;
    else
        data = obj.data.analogData(idx,:);
    end
    sc = 1000;
end
if(~isempty(Args.RemoveLineNoise))
        data = nptRemoveLineNoise(data,Args.RemoveLineNoise,sRate);
end

% set axis
y0 = min(data);
y1 = max(data);

df = length(idx);

if (~isempty(gca))
    delete(gca)
end

hold on
    
for i = 1:length(swrInfo)
    if swrInfo(i,3) >= idx(1) && swrInfo(i,2) <= idx(end) % detects if swr event is in the graph
%             si = 0; %swr event
        if swrInfo(i,3) >= idx(1) && swrInfo(i,2) < idx(1)
            x0 = idx(1);
            x1 = swrInfo(i,3)+0.5;
%                 if swrInfo(i,1) >= idx(1)
%                     si = swrInfo(i,1);
%                 end
        elseif swrInfo(i,2) <= idx(end) && swrInfo(i,3) > idx(end)
            x0 = swrInfo(i,2)-0.5;
            x1 = idx(end);
%                 if swrInfo(i,1) <= idx(end)
%                     si = swrInfo(i,1);
%                 end
        else
            x0 = swrInfo(i,2)-0.5;
            x1 = swrInfo(i,3)+0.5;
%                 si = swrInfo(i,1);
        end
%             si = (si-idx(1)).*(1000/sc);
        patch(([x0 x1 x1 x0]-idx(1)).*(1000/sc), [y0*[1 1] y1*[1 1]], [0.7 0.8 1], 'LineStyle','none')
%             if (si ~= 0)
%                 plot(si, data(si/(1000/sc)),'rx')
%             end
    end
end

plot((1:df)*(1000/sc), data,'Color',Args.Color);

% plotting std away from mean
plot(xlim,[swrMean+3*swrStd swrMean+3*swrStd],'r');
plot(xlim,[swrMean+swrStd swrMean+swrStd],'r'); 
text(1*(1000/sc)+5, 1+swrMean+3*swrStd, 'threshold = 3\sigma');
text(1*(1000/sc)+5, 1+swrMean+swrStd, '1\sigma');
p = gco;
ylim([y0 y1])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% need to move this in due to key input
% @vmswr/PLOT takes 'LabelsOff'
if(~Args.LabelsOff)
    xlabel('time (ms)')
    ylabel('LFP')
end

if(~Args.TitleOff)
    [a,b] = fileparts(obj.nptdata.SessionDirs{:});
    title(b)
    sp = {' '}; %space
    title(strcat('LFP of',sp,b,sp,'from',sp,num2str(idx(1)),'ms to',sp,num2str(idx(end)),'ms'))
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
