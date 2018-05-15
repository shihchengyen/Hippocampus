function placeselect
	
nShuffles = 1000;
res = 5; % 5 x 5 grid resolution
overalGridSize = 25;
oGS2 = overalGridSize/2;
sampleRate = 30000;

numGrids = res*res;
gridSize = overalGridSize/res;
horGridBound = zeros(1,res+1);
vertGridBound = horGridBound;

for e = 1:res+1, % get grid boundaries
    horGridBound(e) = -oGS2 + (e-1)*gridSize; 
    vertGridBound(e) = oGS2 - (e-1)*gridSize;
end

um = unitymaze('auto');
rp = rplparallel('auto');

% data = arrangeMarkers(rh.data,rp);
% data.numSets = size(data.trialIndices,1);
l = load('hmmsort.mat');

unityTriggers = um.data.unityTriggers;
processTrials = um.data.processTrials;
unityData = um.data.unityData;

trigger = rp.data.markers; 
triggerTimes = rp.data.timeStamps; 
pport = vertcat(trigger(1,1:2:end),triggerTimes(1,1:2:end)); % concatenate triggers and timestamps (remove '0's i.e. even columns)
if pport(1,1) > 80, % delete first trigger (trigger version)
    pport(:,1) = [];
end
%rplTriggers(:,1) = pport(2,1:3:end)'; rplTriggers(:,2) = pport(2,2:3:end)'; rplTriggers(:,3) = pport(2,3:3:end)'; rplTrialType = pport(1,3:3:end)';
rplTriggers(:,1) = pport(2,pport(1,:)>10 & pport(1,:)<20)'; 
rplTriggers(:,2) = pport(2,pport(1,:)>20 & pport(1,:)<30)'; 
rplTriggers(:,3) = pport(2,pport(1,:)>30)'; 
rplTrialType = pport(1,pport(1,:)>30)'; % note: exits loop if there are missing ripple triggers for the session as it causes misalignment in matrix
rplTriggers(:,1) = rplTriggers(:,1)*1000; 
rplTriggers(:,2) = rplTriggers(:,2)*1000; 
rplTriggers(:,3) = rplTriggers(:,3)*1000;
rplTrialDur = rplTriggers(:,3) - rplTriggers(:,2);

if(~isempty(l.mlseq))
    ncells = size(l.mlseq,1);
    for si = 1:ncells
	    f1 = figure('name','spike waveform');
	    subplot(2,7,1)
	    plot(l.spikeForms(si,:)); % plot spike form for inspection
	    title('Waveform');

        [wfmin,wfmidx] = min(l.spikeForms(si,:)); % get index of peak
        spikeIdx = find(l.mlseq(si,:) == wfmidx); 
		spikeTimes = spikeIdx/sampleRate; % get peak times
        % Get maze coordinates for each spike during trial (cue offset to end trial)

        for f = 1:nShuffles, % NUMBER OF SHUFFLES
           	trialcount = 0; % count total number of trials
            spikecountSess = 0; % count number of spikes (per trial)
            
            for g = 1:size(unityTriggers,1),   
                if ismember(g,processTrials), % only include one-hit trials
                    trialcount = trialcount + 1;
                    temp = find(spikeTimes > rplTriggers(g,2)/1000 & spikeTimes < rplTriggers(g,3)/1000); % get spikes indices from cue offset to end of trial for this specific trial
                    trialSpkTimes = spikeTimes(temp); 
					trialSpkTimes = trialSpkTimes - rplTriggers(g,2)/1000; % get spikes times align to cue offset time

                    if f > 1,
                        %Do spike shuffling
                        trialSpkTimes = [zeros(1,1),trialSpkTimes]; % pad zero in front
                        orderedISI = diff(trialSpkTimes); % get inter-spike intervals
                        shuffledISI = orderedISI(1,randperm(size(orderedISI,2))); % shuffle ISI
                        shuffledSpikes = cumsum(shuffledISI);
                    else  % if f > 1,
                        shuffledSpikes = trialSpkTimes;
                    end  % if f > 1,

                        unityTrialData = unityData(unityTriggers(g,2):unityTriggers(g,3),:); % extract unity data from cue offset to end of trial

                    % Get location bin for each unity data frame
                    for h = 1:size(unityTrialData,1),
                        xCoord = find(horGridBound < unityTrialData(h,3)); 
						xCoord = xCoord(end);
                        yCoord = find(vertGridBound > unityTrialData(h,4)); 
						yCoord = yCoord(end);
                        unityTrialData(h,6) = xCoord + (yCoord-1)*res;
                    end  % for h = 1:size(unityTrialData,1),

                    % Get duration in each bin per trial
                    for i = 1:numGrids,
						uTDidx = find(unityTrialData(:,6)==i);
						if(~isempty(uTDidx))
							binSpec(trialcount,i) = sum(unityTrialData(uTDidx,2));
						else
							binSpec(trialcount,i) = 0;
						end
                    end  % for i = 1:numGrids,

                    % Get x- y- maze coordinates for each spike during the trial
                    unityTrialTime = cumsum(unityTrialData(:,2)); % get cumulative time from cue offset to end of trial for this specific trial

                    if abs(unityTrialTime(end)-rplTrialDur(g)/1000)< 2; % debugging, ignore all trials with trigger discrepancies > 2 seconds
                        spikecountTrial = 0; % reset trial spike counters
                        clear spikeLocTrial;
                        for j = 1:size(shuffledSpikes,2),
                            spikecountTrial = spikecountTrial + 1; % grow the matrix such that spike locations across the trial is stored in one matrix for the heatmap 
                            temp2 = find(unityTrialTime > shuffledSpikes(1,j),1);
                            spikeLocTrial(spikecountTrial,1:3) = unityTrialData(temp2,[3 4 6]); % matrix for spike locations on the present trial
                        end  % for j = 1:size(shuffledSpikes,2),
                        
                        for aa = 1:25,
                            spaceInd = find(spikeLocTrial(:,3) == aa);
                            anovaMatrix(trialcount,aa) = size(spaceInd,1)/binSpec(trialcount,aa);
                        end  % for aa = 1:25,
                        
                        if f == 1,
                            anovaMatrix1 = anovaMatrix;
                            spikeFreq(f,1:25) = nanmean(anovaMatrix1,1);
                            
                            subplot(2,7,[4 5 11 12])
                            plot(unityTrialData(:,3),unityTrialData(:,4),'Color',[0.5 0.5 0.5],'LineWidth',0.2); 
                            xlim([-12.5 12.5]); 
                            ylim([-12.5 12.5]); 
                            hold on % plot trajectory overlaid with spikes
                            plot(spikeLocTrial(:,1),spikeLocTrial(:,2),'r.','MarkerSize',5); 
                            % hold on
                            set(gca,'XTick',horGridBound); 
                            set(gca,'YTick',fliplr(vertGridBound));
                            grid on
                            title('Trajectory');
                        else  % if f == 1,
                            spikeFreq(f,1:25) = nanmean(anovaMatrix,1);
                        end  % if f == 1,
                    else  % if abs(unityTrialTime(end)-rplTrialDur(g)/1000)< 2; % debugging, ignore all trials with trigger discrepancies > 2 seconds
                        disp('skipped trial');
                        disp(trialcount);
                    end  % if abs(unityTrialTime(end)-rplTrialDur(g)/1000)< 2; % debugging, ignore all trials with trigger discrepancies > 2 seconds
                    
                    % Get average spike frequency for each target
                    spikeTarget(trialcount,1) = size(trialSpkTimes,2); % number of spikes across entire trial (target analysis)
                    spikeTarget(trialcount,2) = spikeTarget(trialcount,1)/unityTrialTime(end); % spike frequency across trial
                    spikeTarget(trialcount,3) = unityTrialData(end,1);              
                end  % if ismember(g,processTrials), % only include one-hit trials
            end  % for g = 1:size(unityTriggers,1),  
        disp(f);
        end  % for f = 1:nShuffles, % NUMBER OF SHUFFLES

        % Get spiking frequency per target location    
        for p = 31:36,
            spikeTargetFreq(p-30,1) = mean(spikeTarget(spikeTarget(:,3) == p,2)); % mean firing rate for the target
            spikeTargetFreq(p-30,2) = std(spikeTarget(spikeTarget(:,3) == p,2)); % mean firing rate for the target
            spikeTargetFreq(p-30,3) = length(find(spikeTarget(:,3) == p)); % number of times target was presented during the session
        end  % for p = 31:36,

        subplot(2,7,8) % Plot target specificity
        bar(spikeTargetFreq(:,1),'FaceColor','green');
        set(gca,'xticklabel',{'T1','T2','T3','T4','T5','T6'});
        title('Target-specificity');
        ylabel('Firing frequency (Hz)'); hold on
        errorbar(1:6,spikeTargetFreq(:,1),spikeTargetFreq(:,2),'k.');
        for q=1:6
            text(q,1,num2str(spikeTargetFreq(q,3)),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',5);
        end 
        
        % pre-ANOVA processing
        samples = 1;
        removeBin = 1;
        for i = 1:25
            timesEntered = anovaMatrix1(~isnan(anovaMatrix1(:,i)),i);
            if length(timesEntered) > 5, % monkey must have traversed that segment at least 5 times throughout the session
                firingRate(samples:samples+length(timesEntered)-1) = timesEntered;
                binNo(samples:samples+length(timesEntered)-1) = i;
                samples = samples+length(timesEntered);
            else  % if length(timesEntered) > 5, % monkey must have traversed that segment at least 5 times throughout the session
                anovaMatrix1(:,i) = nan;
                excludeBins(removeBin) = i;
                removeBin = removeBin + 1;
            end  % if length(timesEntered) > 5, % monkey must have traversed that segment at least 5 times throughout the session
        end  % for i = 1:25
		
		% change to the session directory to read the spreadsheet
		[p,cwd] = getDataOrder('session','CDNow');
        save('excludeBins','excludeBins');
		% return to previous directory
		cd(cwd);
        spikeFreq(:,excludeBins) = nan; % exclude segments traversed < 6 times from heat map analyses
        
        % Get 95th percentile (significance using shuffling method)
        unit_locsig(1,1:numGrids) = spikeFreq(1,:);
        unit_locsig(2,1:numGrids) = prctile(spikeFreq(2:end,1:numGrids),2.5); 
        unit_locsig(3,1:numGrids) = prctile(spikeFreq(2:end,1:numGrids),97.5);

        % Calculate selectivity index (SI)
        maxFiring = max(spikeFreq(1,:));
        fsum = 0;
        for m = 1:25,
            if isnan(spikeFreq(1,m)),
                continue
            end
            fsum = fsum + spikeFreq(1,m)/maxFiring;
        end  % for m = 1:25,
        unit_locSI = (25-length(excludeBins)-fsum)/(25-length(excludeBins));
        disp(unit_locSI);

        % Plot firing rate map
        subplot(2,7,[6 7 13 14])
        firRate = reshape(spikeFreq(1,:),[res,res])'; % reshape spikeFreq matrix to map onto physical maze locations
        imAlpha=ones(size(firRate)); % plot NaN regions as black
        imAlpha(isnan(firRate))=0;
        imagesc(firRate,'AlphaData',imAlpha); % input argument to set color range limits e.g. [100 200]
        set(gca,'color',[0 0 0],'xtick',[],'ytick',[]); % set background to black and remove x- and y-axis labels
        colorbar; hold on % show colorbar
        plot(1.5,2,'r.','MarkerSize',40); % plot poster positions
        plot(4,1.5,'r.','MarkerSize',40);
        plot(4,2.5,'r.','MarkerSize',40);
        plot(2,3.5,'r.','MarkerSize',40);
        plot(2,4.5,'r.','MarkerSize',40);
        plot(4.5,4,'r.','MarkerSize',40); hold on

        for n = 1:25, % highlight significant grids from shuffle
        	row = floor(n/5); 
        	column = rem(n,5); 
            if column == 0, % if no remainder, plot in 5th grid
                column = 5;
            else  % if column == 0, % if no remainder, plot in 5th grid
                row = row + 1;
            end  % if column == 0, % if no remainder, plot in 5th grid
                
            if unit_locsig(1,n) < unit_locsig(2,n), % below cutoff
                textHandles(column,row) = text(column,row,'*','Color','k','FontSize',15,'horizontalAlignment','center'); hold on % indicate significant grids with asterisk
            elseif unit_locsig(1,n) > unit_locsig(3,n), % above cutoff
                textHandles(column,row) = text(column,row,'*','Color','w','FontSize',15,'horizontalAlignment','center'); hold on % indicate significant grids with asterisk
            end  % if unit_locsig(1,n) < unit_locsig(2,n), % below cutoff
        end  % for n = 1:25, % highlight significant grids from shuffle
        
        title('Firing rate map (spikes/s)');
        set(gcf, 'InvertHardCopy', 'off'); % retain black background in figure
        h = suptitle(strcat(getDataOrder('ShortName'), ' Unit ', num2str(si,'%02d'), ' SI: ', num2str(round(unit_locSI,2)))); % main title
        set(h,'FontSize',16);

        set(gcf,'units','normalized','outerposition',[0 0 1 1]); % maximize figure
        set(gcf, 'PaperPositionMode', 'auto');

        % ANOVA and multiple comparisons
        figure();
        [P,T,STATS,TERMS]=anovan(firingRate,{binNo},'display','off'); % use anovan for unbalanced groups
        [COMPARISON,MEANS,H,GNAMES] = multcompare(STATS,'CType','bonferroni','display','off'); % post-hoc pairwise comparisons
        
        figure(f1);
        anovaSig = COMPARISON(COMPARISON(:,6)<0.05,1:2); % find pairwise comparison rows with p < .05
        anovaSigLocs = GNAMES(anovaSig); % get location bins with p < .05
        for bb = 1:numel(anovaSigLocs),
            anovaSigLocs{bb} = anovaSigLocs{bb}(4:end); % remove first 3 characters (redundant for display)
        end  % for bb = 1:numel(anovaSigLocs),
        anovaSigLocs = strcat(anovaSigLocs(:,1),', ',anovaSigLocs(:,2));          
        text(9,3,{'One-way ANOVA'; ' '; strcat('p-value: ',num2str(T{26})); 'post-hoc comparisons: '; char(anovaSigLocs)},'FontSize',8); % display anova statistics on plot

		ustr = ['unit' num2str(si,'%02d') '_placefield'];
	    saveas(gcf,[ustr '.png']);
		save([ustr '.mat'],'unit_locSI')
	    close all;
    
	    clearvars -except nShuffles res numGrids gridSize horGridBound vertGridBound l sampleRate unityTriggers processTrials unityData rplTriggers rplTrialType rplTrialDur spikeSort spikeForms ncells 
    end  % for si = 1:ncells
end  % if(~isempty(l.mlseq))
