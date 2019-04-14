function [obj, varargout] = plot(obj,varargin)
%@vmhighpass/plot Plot function for vmhighpass object.
%   OBJ = plot(OBJ) creates a raster plot of the neuronal
%   response.

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
			'PreTrial',500, 'NormalizeTrial',0, 'RewardMarker',3, ...
            'TimeOutMarker',4, 'PlotAllData',0, 'SpikeData',[], ...
            'SpikeTriggerIndex', 26, 'SpikeHeight', 100, ...
            'FreqPlot',0, 'RemoveLineNoise',[], 'LogPlot',0, ...
		    'FreqLims',[], 'TFfft',0, 'TFfftWindow',200, 'TFfftOverlap',150, ...
		    'TFfftPoints',256, ...
		    'TFWavelets',0,  ...
		    'LoadSort',0, ...
		    'ReturnVars',{''}, 'ArgsOnly',0);
Args.flags = {'LabelsOff','ArgsOnly','NormalizeTrial','FreqPlot','TFfft', ...
				'LogPlot','TFWavelet','PlotAllData','LoadSort'};
[Args,varargin2] = getOptArgs(varargin,Args);

% if user select 'ArgsOnly', return only Args structure for an empty object
if Args.ArgsOnly
    Args = rmfield (Args, 'ArgsOnly');
    varargout{1} = {'Args',Args};
    return;
end

% analogTime needs to be in double-precision in order to avoid rounding
% error, but we don't really want to save it on disk as it takes up too
% much space and makes loading the object slower. So we will create it the
% first time we need it
% check if the field analogTime exists OR if the analogTime field exists
% and is not a double precision variable
% we will create the variable here.
if( ~isfield(obj.data,'analogTime') || ~isa(obj.data.analogTime,'double') )
	obj.data.analogTime = (0:(obj.data.analogInfo.NumberSamples-1))' ./ obj.data.analogInfo.SampleRate;
end

if(~isempty(Args.NumericArguments))
	% plot one data set at a time
	n = Args.NumericArguments{1};
	% find indices for n-th trial
	tIdx = obj.data.trialIndices(n,:);
	if(length(tIdx)==2)
		OldMarkerFormat = 1;
	else
		OldMarkerFormat = 0;
	end
	sRate = obj.data.analogInfo.SampleRate;
	if(OldMarkerFormat)
		idx = (tIdx(1)-(Args.PreTrial/1000*sRate)):tIdx(2);
	else
		idx = (tIdx(1)-(Args.PreTrial/1000*sRate)):tIdx(3);
	end
	if(Args.FreqPlot)
		if(Args.PlotAllData)
			data = obj.data.analogData;
		else
			data = obj.data.analogData(idx);
		end
		if(~isempty(Args.RemoveLineNoise))
			data = nptRemoveLineNoise(data,Args.RemoveLineNoise,sRate);
		end
		% remove the mean, i.e. DC component		
		datam = mean(data);			
		PlotFFT(data-datam,sRate);
		set(gca,'TickDir','out')
		if(Args.LogPlot)
			set(gca,'YScale','log')
		end
	elseif(Args.TFfft)
		if(Args.PlotAllData)
			% create memory to store overall mean
			if(OldMarkerFormat)
				dIdx = diff(obj.data.trialIndices,1,2);
			else
				dIdx = obj.data.trialIndices(:,3) - obj.data.trialIndices(:,1);
			end
			% find longest trial
			mIdx = max(dIdx);
			% find number of time bins in the spectrogram that corresponds to
			spTimeStep = Args.TFfftWindow - Args.TFfftOverlap;
			spTimeBins = floor(mIdx/spTimeStep) - Args.TFfftOverlap/spTimeStep;
			% create matrix
			nFreqs = (Args.TFfftPoints/2)+1;			
			ops = zeros(nFreqs,spTimeBins);
			opsCount = ops;
			for ti = 1:obj.data.numSets
				tftIdx = obj.data.trialIndices(ti,:);
				if(OldMarkerFormat)
					tfidx = tftIdx(1):tftIdx(2);
				else
					tfidx = tftIdx(1):tftIdx(3);
				end
				data = obj.data.analogData(tfidx);
				if(~isempty(Args.RemoveLineNoise))
					data = nptRemoveLineNoise(data,Args.RemoveLineNoise,sRate);
				end
				datam = mean(data);
				[s,w,t,ps] = spectrogram(data-datam,Args.TFfftWindow, ...
					Args.TFfftOverlap,Args.TFfftPoints,sRate);
				% add to overall mean
				% get columns of ps
				psIdx = 1:size(ps,2);
				ops(:,psIdx) = ops(:,psIdx) + ps;
				opsCount(:,psIdx) = opsCount(:,psIdx) + 1;
			end
			imagesc(0:(Args.TFfftWindow-Args.TFfftOverlap):mIdx,0:(sRate/Args.TFfftPoints):(sRate/2),ops./opsCount)
            set(gca,'Ydir','normal')
		else
			data = obj.data.analogData(idx);
			if(~isempty(Args.RemoveLineNoise))
				data = nptRemoveLineNoise(data,Args.RemoveLineNoise,sRate);
			end
			datam = mean(data);
			spectrogram(data-datam,Args.TFfftWindow,Args.TFfftOverlap,Args.TFfftPoints, ...
				sRate,'yaxis')
		end
	elseif(Args.TFWavelets)
		% not fully completed yet
    	cdata = nptRemoveLineNoise(obj.data.analogData',50,1000);
    	cdata1 = cdata(idx);
    	cdata1m = mean(cdata1);
    	PlotFFT(cdata1-cdata1m,1000);
    
        data.trial = num2cell(mdata,2);
        data.label = '1';
        data.fsample = 1000;
        data.time = num2cell(repmat((0:13087)/1000,51,1),2);

        cfg.channel = 'all';
        cfg.method = 'wavelet';
        cfg.width = 7;
        cfg.output = 'pow';
        cfg.foi = 1:100;
        cfg.toi = 0:5000;
        cfg.pad = 'nextpow2';
        
        TFRwave = ft_freqanalysis(cfg,data);
	else
		data = obj.data.analogData(idx);
		if(~isempty(Args.RemoveLineNoise))
			data = nptRemoveLineNoise(data,Args.RemoveLineNoise,sRate);
		end
		plot( (obj.data.analogTime(idx)-obj.data.analogTime(tIdx(1))) * 1000,data,'.-')
		% indicate trial start
		line([0 0],ylim,'Color','g')
		if(OldMarkerFormat)
			if(obj.data.markers(n,2)==Args.RewardMarker)
				% indicate correct trial
				line(repmat((obj.data.analogTime(idx(end))-obj.data.analogTime(tIdx(1)))*1000,2,1),ylim,'Color','b')
			else
				% indicate incorrect trial
				line(repmat((obj.data.analogTime(idx(end))-obj.data.analogTime(tIdx(1)))*1000,2,1),ylim,'Color','r')
			end
		else
			% indicate cue offset
			line(repmat((obj.data.analogTime(tIdx(2))-obj.data.analogTime(tIdx(1)))*1000,2,1),ylim,'Color','m')
			if(obj.data.markers(n,3)==Args.RewardMarker || (floor(obj.data.markers(n,3)/10)==Args.RewardMarker) )
				% indicate correct trial
				line(repmat((obj.data.analogTime(idx(end))-obj.data.analogTime(tIdx(1)))*1000,2,1),ylim,'Color','b')
			else
				% indicate incorrect trial
				line(repmat((obj.data.analogTime(idx(end))-obj.data.analogTime(tIdx(1)))*1000,2,1),ylim,'Color','r')
			end
		end
		if(Args.LoadSort)
            % if UserData does not contain hmmsort, then save it there so
            % we don't have to load it every trial
            a = get(gcf,'UserData');
            if(~isfield(a,'hmmsort'))
                l = load('hmmsort.mat');
                a.hmmsort = l;
                set(gcf,'UserData',a);
            else
                l = a.hmmsort;
            end
            Args.SpikeData = l.mlseq;
		end			
		if(~isempty(Args.SpikeData))
			hold on
			% get SpikeData
			mlseq = Args.SpikeData;
			spidx = Args.SpikeTriggerIndex;
			% add spike trains
            ncells = size(mlseq,1);
       		clist = nptDefaultColors(1:ncells);
            % first get size of data being plotted
            idxsize = size(idx,2);
            % create a vector that is the same size
            a = nan(idxsize,1);
            % create a vector that will be used to plot overlapping waveforms
            as = zeros(idxsize,1);
            % keep track of which indices need to be plotted because they
            % are overlaps
            asi = as;
            % get size of waveforms
            wfsize = size(l.spikeForms,3);
            % create vector from 1 to wfsize
            wfsvec = 1:wfsize;
            for si = 1:ncells
                st1 = find(mlseq(si,idx)==spidx);
                if(~isempty(st1))
                    % add stem plot
                    stem( (obj.data.analogTime(idx(st1))-obj.data.analogTime(tIdx(1))) * 1000, repmat(Args.SpikeHeight,[size(st1),1]), 'Color', clist(si,:))
                    % plot the templates
                    % find the starting index for each spike
                    astart = st1-spidx;
                    % get number of spikes
                    nspikes = size(st1,2);
                    % create array of indices
                    aind = repmat(wfsvec,nspikes,1)' + repmat(astart,wfsize,1);
                    % check to make sure aind does not contain indices 
                    % larger than idxsize or smaller than 1
                    aind(aind>idxsize) = idxsize;
                    aind(aind<1) = 1;
                    % create spike waveforms
                    swaveforms = repmat(squeeze(l.spikeForms(si,:,:)),1,nspikes);
                    % initialize ai to all nan's
                    ai = a;
                    ai(aind) = swaveforms;
		            plot( (obj.data.analogTime(idx)-obj.data.analogTime(tIdx(1))) * 1000, ai,'Color', clist(si,:), 'LineStyle',':');
                    % insert the spike waveform, and add to existing values
                    % for overlapping waveforms
                    as(aind) = as(aind) + swaveforms;
                    asi(aind) = asi(aind) + 1;
                end  % if(~isempty(st1))
            end  % for si = 1:ncells
            % find points with overlaps
            asio = asi>1;
            % set those points to 2 since we just need to differentiate
            % them from points with no overlaps
            asi(asio) = 2;
            % find points where we transition from no overlap to overlap
            % and vice versa
            dasi = diff(asi);
            % find points where dasi is positive
            dasip = dasi>0;
            % find points where dasi is negative
            dasim = dasi<0;
            % we need to shift the indices to account for the fact that
            % diff will produce vectors shorter than the original vectors
            dasip2 = [dasip; 0];
            dasim2 = [0; dasim];
            % combine both together
            asit = dasip2 | dasim2;
            % combine with indices with overlaps
            asc = asio | asit;
            % remove the other points
            as(~asc) = nan;
            % plot spike waveforms with dotted lines
            plot( (obj.data.analogTime(idx)-obj.data.analogTime(tIdx(1))) * 1000, as,'Color', 'k', 'LineStyle','--');
			hold off
		end	
	end
else
	% plot all data
	if(Args.FreqPlot)
		sRate = obj.data.analogInfo.SampleRate;
		data = obj.data.analogData;
		if(~isempty(Args.RemoveLineNoise))
			data = nptRemoveLineNoise(data,Args.RemoveLineNoise,sRate);
		end

		datam = mean(data);
		PlotFFT(data-datam,sRate)
		set(gca,'TickDir','out')
		if(Args.LogPlot)
			set(gca,'YScale','log')
		end
	else
		if(size(obj.data.trialIndices,2)<3)
			OldMarkerFormat = 1;
		else
			OldMarkerFormat = 0;
		end
		
		if(OldMarkerFormat)
			dIdx = diff(obj.data.trialIndices,1,2);
		else
			dIdx = obj.data.trialIndices(:,3) - obj.data.trialIndices(:,1);
		end
		% find longest trial
		mIdx = max(dIdx);
		% create matrix
		mdata = zeros(obj.data.numSets,mIdx);
		for i = 1:obj.data.numSets
			idx = obj.data.trialIndices(i,:);
			if(Args.NormalizeTrial)
				if(ldMarkerFormat)
					rdata = obj.data.analogData(idx(1):idx(2));
				else
					rdata = obj.data.analogData(idx(1):idx(3));
				end
				rdmin = min(rdata);
				rdmax = max(rdata);
				mdata(i,1:(dIdx(i)+1)) = (rdata-rdmin)/(rdmax-rdmin);
			else
				if(OldMarkerFormat)
					mdata(i,1:(dIdx(i)+1)) = obj.data.analogData(idx(1):idx(2));
				else
					mdata(i,1:(dIdx(i)+1)) = obj.data.analogData(idx(1):idx(3));
				end
			end
		end
		imagesc(mdata)
		colormap(jet)
	end
end

% add code for plot options here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% @vmhighpass/PLOT takes 'LabelsOff' as an example
if(~Args.LabelsOff)
	if(Args.FreqPlot)
		xlabel('Frequency (Hz)')
		ylabel('Magnitude')
	elseif(Args.TFfft)
		xlabel('Time (s)')
		ylabel('Frequency (Hz)')
	else
		xlabel('Time (ms)')
		ylabel('Voltage (uV)')
	end
end
sdstr = get(obj,'SessionDirs');
title(getDataOrder('ShortName','DirString',sdstr{1}))
if(~isempty(Args.FreqLims))
	if(Args.FreqPlot)
		xlim(Args.FreqLims)
	elseif(Args.TFfft)
		ylim(Args.FreqLims)
	end
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
