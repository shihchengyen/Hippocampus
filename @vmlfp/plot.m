function [obj, varargout] = plot(obj,varargin)
%@vmlfp/plot Plot function for vmlfp object.
%   OBJ = plot(OBJ) creates a raster plot of the neuronal
%   response.

Args = struct('LabelsOff',0,'GroupPlots',1,'GroupPlotIndex',1,'Color','b', ...
			'PreTrial',500, 'NormalizeTrial',0, 'RewardMarker',3, ...
            'TimeOutMarker',4, 'PlotAllData',0, 'TitleOff', 0, ...
            'FreqPlot',0, 'RemoveLineNoise',[], 'LogPlot',0, ...
		    'FreqLims',[], 'TFfft',0, 'TFfftWindow',200, 'TFfftOverlap',150, ...
		    'TFfftPoints',256, 'TFfftStart',500,'TFfftFreq',150,...
		    'TFWavelets',0, 'TimeWindow', [],  ...
            'Filter',0,'OverlapLFP',0,'FilterWindow',[],'CorrCoeff',0,   ...
		    'ReturnVars',{''}, 'ArgsOnly',0);
Args.flags = {'LabelsOff','ArgsOnly','NormalizeTrial','FreqPlot','TFfft', ...
				'LogPlot','TFWavelet','PlotAllData','TitleOff','Filter','OverlapLFP','CorrCoeff'};
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
    
    % =====================================================================
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
        
    % =====================================================================    
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
            tIdx = obj.data.trialIndices(n,:); % Obtain trial indexes
            
            %     Spectrogram data for the 'normalisation period NP'
            idx = (tIdx(1)-((Args.TFfftStart+500)/1000*sRate)):(tIdx(1)-((Args.TFfftStart+1)/1000*sRate));
            %Inter-trial interval data used for normalisation will always
            %be 500ms before TfftStart
            data = obj.data.analogData(idx);
            datam = mean(data);
            [~,~,~,P]=spectrogram(data-datam,Args.TFfftWindow,Args.TFfftOverlap,Args.TFfftPoints,sRate,'yaxis');
            
            %     Normalization parameters of the NP
            Pmean=mean(P,2); %mean power density of each frequency bin
            Pstd=std(P,0,2); %standard deviation of each frequency bin
            
            %     Spectrogram data for trials
            idx = (tIdx(1)-(Args.TFfftStart/1000*sRate)):tIdx(3);
            %Trial data including pre-trial data of duration determined by TfftStart
            data = obj.data.analogData(idx);
            datam = mean(data);
            [spec.S,spec.F,spec.T,spec.P,spec.Fc,spec.Tc]=...
                spectrogram(data-datam,Args.TFfftWindow,Args.TFfftOverlap,Args.TFfftPoints,sRate,'yaxis');
            
            % trial psd normalised to NP
            spec.Pnorm=(spec.P-Pmean)./Pstd;
            spec.T=(-Args.TFfftStart/1000:(Args.TFfftWindow-Args.TFfftOverlap)/sRate:spec.T(end)-(Args.TFfftStart/1000+Args.TFfftWindow/sRate/2));
           
            % Plots the normalized PSD data in a figure
            surf(spec.T,spec.F,spec.Pnorm,'EdgeColor','none');
            axis xy; axis([-Args.TFfftStart/1000 inf 0 Args.TFfftFreq]); colormap(jet); view(0,90); caxis([-10 10]);
            set(gca,'FontSize',6); xticks(-Args.TFfftStart/1000:0.5:spec.T(end)); yticks(0:10:Args.TFfftFreq);
            title(strcat('Normalised Spectrogram of Trial:',string(n)),'FontSize',10);
            colorbar;
            
            % Plot lines to mark the cue presentation period
            hold on
            line([0 0],[0 150],[100 100; 100 100],'Color','k');
            line([1 1],[0 150],[100 100; 100 100],'Color','k');
            hold off
            
% 			data = obj.data.analogData(idx);
% 			if(~isempty(Args.RemoveLineNoise))
% 				data = nptRemoveLineNoise(data,Args.RemoveLineNoise,sRate);
% 			end
% 			datam = mean(data);
% 			spectrogram(data-datam,Args.TFfftWindow,Args.TFfftOverlap,Args.TFfftPoints, ...
% 				sRate,'yaxis','MinThreshold',-10)
%             ylim([0,100])
        end
    
    % =====================================================================    
    elseif(Args.CorrCoeff)
        % Computes the correlation coefficient plot with the same
        % normalisation period (-1s to -0.5s) but with NO WINDOW OVERLAPS
        
        tIdx = obj.data.trialIndices(n,:); % Obtain trial indexes
            
            %     Spectrogram data for the 'normalisation period NP'
            idx = (tIdx(1)-((Args.TFfftStart+500)/1000*sRate)):(tIdx(1)-((Args.TFfftStart+1)/1000*sRate));
            data = obj.data.analogData(idx);
            datam = mean(data);
            [~,~,~,P]=spectrogram(data-datam,Args.TFfftWindow,0,Args.TFfftPoints,sRate,'yaxis');
            
            %     Normalization parameters of the NP
            Pmean=mean(P,2); %mean power density of each frequency bin
            Pstd=std(P,0,2); %standard deviation of each frequency bin
            
            %     Spectrogram data for trials
            idx = (tIdx(1)-(Args.TFfftStart/1000*1000)):tIdx(3);
            data = obj.data.analogData(idx);
            datam = mean(data);
            [spec.S,spec.F,spec.T,spec.P,spec.Fc,spec.Tc]=...
                spectrogram(data-datam,Args.TFfftWindow,0,Args.TFfftPoints,sRate,'yaxis');
            
            % trial psd normalised to NP
            spec.Pnorm=(spec.P-Pmean)./Pstd;

            %     Computes correlation coefficients for each trial
            CCmean=mean(spec.Pnorm,2); %mean power density of each frequency spectrum of trials
            Pdiff=zeros((size(spec.Pnorm,1)),(size(spec.Pnorm,2))); % Initialises Pdiff
            
            for r=1:size(spec.Pnorm,1)
                % Difference at each time point with the mean power density for that frequency
                Pdiff(r,:)=spec.Pnorm(r,:)-CCmean(r);
            end
            
            for r=1:size(spec.Pnorm,1)
                for c=1:size(spec.Pnorm,1)
                    %Correlation coefficient for each frequency band
                    spec.cc(r,c)=sum(Pdiff(r,:).*Pdiff(c,:))/...
                        sqrt(sum(Pdiff(r,:).*Pdiff(r,:),2)*sum(Pdiff(c,:).*Pdiff(c,:),2));
                end
            end
            
            % Determines the data to enter into XLimits and YLimits based
            % on the TFfftFreq limit initialised as an option
            [~,FreqI]=min(abs(spec.F-Args.TFfftFreq));
            FreqLim=spec.F(FreqI);

            % Plots the CC plots in one figure
            heatmap(spec.F,flipud(spec.F),flipud(spec.cc),...
                'XLimits',{0,FreqLim},'YLimits',{FreqLim,0},...
                'XDisplayLabels',(round(spec.F)),'YDisplayLabels',(flipud(round(spec.F))),...
                'XLabel','Frequency (Hz)','YLabel','Frequency (Hz)','FontSize',6,...
                'Colormap', jet,'ColorLimits',[-0.2,0.7],'ColorbarVisible','on',...
                'Title',strcat('Correlation coefficients of Trial:',string(n)));
            
	% =====================================================================
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
        
    % =====================================================================
    elseif (Args.Filter) 
        % Filter is a flag that plots the time-domain signal with a
        % bandpass filter defined by 'FilterWindow'
        data = obj.data.analogData(idx,:);
        drows = size(data,2);
		if(~isempty(Args.RemoveLineNoise))
			data = nptRemoveLineNoise(data,Args.RemoveLineNoise,sRate);
		end
		axesPositions = separateAxis('Vertical',drows);
		for spi = 1:drows
			% do subplot if there are multiple rows of data
			subplot('Position',axesPositions(spi,:));
			plot( (obj.data.analogTime(idx)-obj.data.analogTime(tIdx(1)) )*1000,...
                bandpass(data(:,spi),Args.FilterWindow,sRate),'.-')
            % bandpass function used to filter the time domain signal
            if(spi>1)
                set(gca,'XTickLabel',{})
            end
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
				if( (obj.data.markers(n,3)==Args.RewardMarker) || (floor(obj.data.markers(n,3)/10)==Args.RewardMarker) )
					% indicate correct trial
					line(repmat((obj.data.analogTime(idx(end))-obj.data.analogTime(tIdx(1)))*1000,2,1),ylim,'Color','b')
				else
					% indicate incorrect trial
					line(repmat((obj.data.analogTime(idx(end))-obj.data.analogTime(tIdx(1)))*1000,2,1),ylim,'Color','r')
				end
			end
        end
    
    % ===================================================================== 
    elseif (Args.OverlapLFP) 
        % Filter is a flag that plots the time-domain signal with a
        % bandpass filter defined by 'FilterWindow'
        
        % ==================================================
        % Look at only the first 1s of the navigation period
        idx = (tIdx(1)+(1000/1000*sRate)):(tIdx(1)+(1999/1000*sRate));
        % ==================================================
        
        data = obj.data.analogData(idx,:);
        drows = size(data,2);
		if(~isempty(Args.RemoveLineNoise))
			data = nptRemoveLineNoise(data,Args.RemoveLineNoise,sRate);
		end
		axesPositions = separateAxis('Vertical',drows);
        for spi = 1:drows
            % do subplot if there are multiple rows of data
            subplot('Position',axesPositions(spi,:));
            
        % ==================================================
        % Look at only the first 1s of the navigation period
            plot( (obj.data.analogTime(idx)-obj.data.analogTime(tIdx(1)) )*1000,data(:,spi),'.-')
            hold on
            plot( (obj.data.analogTime(idx)-obj.data.analogTime(tIdx(1)) )*1000,...
                bandpass(data(:,spi),Args.FilterWindow,sRate),'.-','Color','r')
            hold off
        % ==================================================
        
%             plot( (obj.data.analogTime(idx)-obj.data.analogTime(tIdx(1)) )*1000,data(:,spi),'.-')
%             hold on
%             plot( (obj.data.analogTime(idx)-obj.data.analogTime(tIdx(1)) )*1000,...
%                 bandpass(data(:,spi),Args.FilterWindow,sRate),'.-','Color','r')
%             hold off
           
            % bandpass function used to filter the time domain signal
            
            if(spi>1)
                set(gca,'XTickLabel',{})
            end
            % indicate trial start
%             line([0 0],ylim,'Color','g')
%             if(OldMarkerFormat)
%                 if(obj.data.markers(n,2)==Args.RewardMarker)
%                     % indicate correct trial
%                     line(repmat((obj.data.analogTime(idx(end))-obj.data.analogTime(tIdx(1)))*1000,2,1),ylim,'Color','b')
%                 else
%                     % indicate incorrect trial
%                     line(repmat((obj.data.analogTime(idx(end))-obj.data.analogTime(tIdx(1)))*1000,2,1),ylim,'Color','r')
%                 end
%             else
%                 % indicate cue offset
%                 line(repmat((obj.data.analogTime(tIdx(2))-obj.data.analogTime(tIdx(1)))*1000,2,1),ylim,'Color','m')
%                 if( (obj.data.markers(n,3)==Args.RewardMarker) || (floor(obj.data.markers(n,3)/10)==Args.RewardMarker) )
%                     % indicate correct trial
%                     line(repmat((obj.data.analogTime(idx(end))-obj.data.analogTime(tIdx(1)))*1000,2,1),ylim,'Color','b')
%                 else
%                     % indicate incorrect trial
%                     line(repmat((obj.data.analogTime(idx(end))-obj.data.analogTime(tIdx(1)))*1000,2,1),ylim,'Color','r')
%                 end
%             end
        end
        
    % =====================================================================    
    else
		data = obj.data.analogData(idx,:);
        drows = size(data,2);
		if(~isempty(Args.RemoveLineNoise))
			data = nptRemoveLineNoise(data,Args.RemoveLineNoise,sRate);
		end
		axesPositions = separateAxis('Vertical',drows);
		for spi = 1:drows
			% do subplot if there are multiple rows of data
			subplot('Position',axesPositions(spi,:));
			plot( (obj.data.analogTime(idx)-obj.data.analogTime(tIdx(1)) )*1000,data(:,spi),'-')
            if(spi>1)
                set(gca,'XTickLabel',{})
            end
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
				if( (obj.data.markers(n,3)==Args.RewardMarker) || (floor(obj.data.markers(n,3)/10)==Args.RewardMarker) )
					% indicate correct trial
					line(repmat((obj.data.analogTime(idx(end))-obj.data.analogTime(tIdx(1)))*1000,2,1),ylim,'Color','b')
				else
					% indicate incorrect trial
					line(repmat((obj.data.analogTime(idx(end))-obj.data.analogTime(tIdx(1)))*1000,2,1),ylim,'Color','r')
				end
			end
			
			if(~isempty(Args.TimeWindow))
				xlim(Args.TimeWindow)
			end
		end
    end

% ========================================================================    
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
				if(OldMarkerFormat)
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
% @vmlfp/PLOT takes 'LabelsOff' as an example
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
if(~Args.TitleOff)
    [a,b] = fileparts(obj.nptdata.SessionDirs{:});
    title(b)
end
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
