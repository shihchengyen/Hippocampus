function [obj, varargout] = plot(obj,varargin)
%@vmlfp/plot Plot function for vmlfp object.
%   OBJ = plot(OBJ) creates a raster plot of the neuronal
%   response.

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
			'PreTrial',500, 'NormalizeTrial',0, 'RewardMarker', 3, 'TimeOutMarker', 4, ...
		  'ReturnVars',{''}, 'ArgsOnly',0);
Args.flags = {'LabelsOff','ArgsOnly','NormalizeTrial'};
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
	% find indices for n-th trial
	tIdx = obj.data.trialIndices(n,:);
	idx = (tIdx(1)-(Args.PreTrial/1000*obj.data.analogInfo.SampleRate)):tIdx(2);
	plot((obj.data.analogTime(idx)-obj.data.analogTime(tIdx(1)))*1000,obj.data.analogData(idx),'.-')
    line([0 0],ylim,'Color','g')
    if(obj.data.markers(n,2)==Args.RewardMarker)
        line(repmat((obj.data.analogTime(idx(end))-obj.data.analogTime(tIdx(1)))*1000,2,1),ylim,'Color','b')
    else
        line(repmat((obj.data.analogTime(idx(end))-obj.data.analogTime(tIdx(1)))*1000,2,1),ylim,'Color','r')
    end
else
	% plot all data
	dIdx = diff(obj.data.trialIndices,1,2);
	% find longest trial
	mIdx = max(dIdx);
	% create matrix
	mdata = zeros(obj.data.numSets,mIdx);
	for i = 1:obj.data.numSets
		idx = obj.data.trialIndices(i,:);
        if(Args.NormalizeTrial)
            rdata = obj.data.analogData(idx(1):idx(2));
            rdmin = min(rdata);
            rdmax = max(rdata);
            mdata(i,1:(dIdx(i)+1)) = (rdata-rdmin)/(rdmax-rdmin);
        else
            mdata(i,1:(dIdx(i)+1)) = obj.data.analogData(idx(1):idx(2));
        end
	end
	imagesc(mdata)
    colormap(jet)
end

% add code for plot options here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% @vmlfp/PLOT takes 'LabelsOff' as an example
if(~Args.LabelsOff)
	xlabel('Time (ms)')
	ylabel('Voltage (uV)')
end
[a,b] = fileparts(obj.nptdata.SessionDirs{:});
title(b)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RR = eval('Args.ReturnVars');
for i=1:length(RR) RR1{i}=eval(RR{i}); end 
varargout = getReturnVal(Args.ReturnVars, RR1);
