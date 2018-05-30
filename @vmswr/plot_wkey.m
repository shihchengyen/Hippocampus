% takes in key input
% the key input does not seem compatible with the current structure

function [obj, varargout] = plot(obj,varargin)
%@vmswr/plot Plot function for vmswr object.
%   OBJ = plot(OBJ) creates a raster plot of the neuronal
%   response.

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
            'PreTrial',500,'Rms',1, ...
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
%     % find indices for n-th trial
% 	tIdx = obj.data.trialIndices(n,:);
%     sRate = obj.data.analogInfo.SampleRate;
%     swrInfo = obj.data.analogInfo.Swr;
%     
% 	if(length(tIdx)==2)
% 		OldMarkerFormat = 1;
% 	else
% 		OldMarkerFormat = 0;
%     end
%     
%     if(OldMarkerFormat)
% 		idx = (tIdx(1)-(Args.PreTrial/1000*sRate)):tIdx(2);
% 	else
% 		idx = (tIdx(1)-(Args.PreTrial/1000*sRate)):tIdx(3);
%     end
%     
%     if(Args.Rms)
%         if (Args.PlotAllData)
%             data = obj.data.analogRmsData;
%         else
%             data = obj.data.analogRmsData(idx,:);
%         end
%         sc = 200; %scaling for graphs
%         numPlots = floor(size(data,2)/sc);
%     else
%         if (Args.PlotAllData)
%             data = obj.data.analogData;
%         else
%             data = obj.data.analogData(idx,:);
%         end
%         sc = 1000;
%         numPlots = floor(size(data,2)/sc);
%     end
%     if(~isempty(Args.RemoveLineNoise))
% 			data = nptRemoveLineNoise(data,Args.RemoveLineNoise,sRate);
%     end
%     
%     j = 1;
%     while j < numPlots+1
%         
%         xi = (j-1)*200+1;
%         xf = j*200;
%         yi = (j-1)*sc+1;
%         yf = j*sc;
%         
%         y0 = min(data(1, yi:yf));
%         y1 = max(data(1, yi:yf));
%         
%         for i = 1:length(swrInfo)
%             if swrInfo(i,3) >= xi || swrInfo(i,2) <= xf % detects if swr event is in the graph
%                 if swrInfo(i,3) >= xi && swrInfo(i,2) < xi
%                     x0 = xi;
%                     x1 = swrInfo(i,3)+0.5;
%                 elseif swrInfo(i,2) <= xf && swrInfo(i,3) > xf
%                     x0 = swrInfo(i,2)-0.5;
%                     x1 = xf;
%                 else
%                     x0 = swrInfo(i,2)-0.5;
%                     x1 = swrInfo(i,3)+0.5;
%                 end
%                 
%                 patch([x0 x1 x1 x0].*5, [y0*[1 1] y1*[1 1]], [0.7 0.8 1], 'LineStyle','none') 
%             end
%         end
%         
%         plot((j-1)*sc+1:j*sc, data(1,(j-1)*sc+1:j*sc)); %highlights the SPW Event
%         hold on
%         plot((j-1)*sc+1:j*sc, data(1,yi:yf));
%     end
    
else
	% plot all data
	n = 1;
    fprintf('work in progress\n');
end

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
    numPlots = floor(size(data,1)/sc);
end
if(~isempty(Args.RemoveLineNoise))
        data = nptRemoveLineNoise(data,Args.RemoveLineNoise,sRate);
end

j = 1; % counter for numPlots
% set axis
y0 = min(data);
y1 = max(data);
    
while j < numPlots+1
    
    if (~isempty(gca))
        delete(gca)
    end
    
    ri = idx(1)+(j-1)*200;
    rf = idx(1)+j*200-1;
    di = (j-1)*sc+1;
    df = j*sc;

%     y0 = min(data(di:df, 1));
%     y1 = max(data(di:df, 1));
    
    hold on
    
    for i = 1:length(swrInfo)
        if swrInfo(i,3) >= ri && swrInfo(i,2) <= rf % detects if swr event is in the graph
%             si = 0; %swr event
            if swrInfo(i,3) >= ri && swrInfo(i,2) < ri
                x0 = ri;
                x1 = swrInfo(i,3)+0.5;
%                 if swrInfo(i,1) >= ri
%                     si = swrInfo(i,1);
%                 end
            elseif swrInfo(i,2) <= rf && swrInfo(i,3) > rf
                x0 = swrInfo(i,2)-0.5;
                x1 = rf;
%                 if swrInfo(i,1) <= rf
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

    plot((di:df)*(1000/sc), data(di:df,1),'Color',Args.Color);
    
    % plotting std away from mean
    plot(xlim,[swrMean+3*swrStd swrMean+3*swrStd],'r');
    plot(xlim,[swrMean+swrStd swrMean+swrStd],'r'); 
    text(di*(1000/sc)+5, 1+swrMean+3*swrStd, 'threshold = 3\sigma');
    text(di*(1000/sc)+5, 1+swrMean+swrStd, '1\sigma');
    p = gco;
    ylim([y0 y1])
    
    % plot label options
    if(~Args.LabelsOff)
        xlabel('time (ms)')
        ylabel('LFP')
    end

    if(~Args.TitleOff)
        [a,b] = fileparts(obj.nptdata.SessionDirs{:});
        title(b)
        sp = {' '}; %space
        title(strcat('LFP of',sp,b,sp,'from',sp,num2str(ri),'ms to',sp,num2str(rf),'ms'))
    end
    
    % key input switch (left or right arrow)
    [~,~,button]=ginput(1);
    switch button
      case 28 %left
          if j > 1
            j = j - 1;
          end
      case 29 %right
          if j < numPlots
              j = j + 1;
          end
    end
    
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % need to move this in due to key input
% % @vmswr/PLOT takes 'LabelsOff'
% if(~Args.LabelsOff)
% 	xlabel('time (ms)')
% 	ylabel('LFP')
% end
% 
% if(~Args.TitleOff)
%     [a,b] = fileparts(obj.nptdata.SessionDirs{:});
%     title(b)
% end
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
