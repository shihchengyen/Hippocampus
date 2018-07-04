function bandselect(dirstr)
    
    %% Load Data
    %Behavioural Data
    um = unitymaze('auto');
    processTrials = um.data.processTrials;
    unityData = um.data.unityData;
    unityTriggers = um.data.unityTriggers;
    
    %Ripple Data
    rp = rplparallel('auto');
    trigger = rp.data.markers;
    triggerTimes = rp.data.timeStamps;
    pport = vertcat(trigger(1,1:2:end),triggerTimes(1,1:2:end)); % concatenate triggers and timestamps (remove 0)
    
    if pport(1,1) > 80 % delete first trigger (trigger version)
        pport(:,1) = [];
    else
    end
    
    rplTriggers(:,1) = pport(2,1:3:end)';
    rplTriggers(:,2) = pport(2,2:3:end)';
    rplTriggers(:,3) = pport(2,3:3:end)';
    rplTriggers(:,1) = rplTriggers(:,1)*1000;
    rplTriggers(:,2) = rplTriggers(:,2)*1000;
    rplTriggers(:,3) = rplTriggers(:,3)*1000;

    %LFP
    rl = load('rpllfp.mat');
    lfpWaveform = rl.df.data.analogData;
    lfpWaveform = nptRemoveLineNoise(lfpWaveform,50,1000); % remove 50Hz line noise
    
    %% Place Cell
    res = 5; % 5 x 5 grid resolution
    numGrids = res*res;
    gridSize = 25/res;
    
    horGridBound = zeros(1,res+1);
    vertGridBound = zeros(1,res+1);

    for n = 1:res+1
        horGridBound(n) = -12.5 + (n-1)*gridSize; % get grid boundaries
        vertGridBound(n) = 12.5 - (n-1)*gridSize;
    end
    
    for p = 1:10 % NUMBER OF SHUFFLES (originally 1000)
        binCounter = 0; % bin counter
        proctrialsCounter = 0; % processed trials counter
        q = 1; % trial counter
        while q <= size(rplTriggers,1) % trial-by-trial
            if ismember(q,processTrials)
                try    
                     lfp1 = lfpWaveform(round(rplTriggers(q+1,2)):round(rplTriggers(q+1,3)));
                     lfp2 = lfpWaveform(round(rplTriggers(q,3)):round(rplTriggers(q+1,1)));
                     [s1,f1,t1,ps1] = spectrogram(lfp1,500,250,512,1000,'yaxis'); % in-trial (start at trial 2 for baseline normalization)
                     [s2,f2,t2,ps2] = spectrogram(lfp2,500,250,512,1000,'yaxis'); % ITI/baseline
                     trialdata = unityData(unityTriggers(q+1,2):unityTriggers(q+1,3),:); trialdata(:,2) = cumsum(trialdata(:,2)); % extract unity data for current trial 

                     if p > 1
                         ps1 = ps1(:,randperm(size(ps1,2))); %Shuffle windows
                         ps2 = ps2(:,randperm(size(ps2,2)));
                     else
                         proctrialsCounter = proctrialsCounter + 1; 
                         % get average power over time during cue and navigation periods
                         theta_nav = mean(ps1(3:5,:)); theta_cue = mean(ps2(3:5,:));
                         alpha_nav = mean(ps1(6:8,:)); alpha_cue = mean(ps2(6:8,:));
                         gamma_nav = mean(ps1(17:62,:)); gamma_cue = mean(ps2(17:62,:));
                         thetaTime(proctrialsCounter,1:size(ps2,2)) = theta_cue; thetaTime(proctrialsCounter,size(ps2,2)+1:size(ps1,2)+7) = theta_nav;
                         alphaTime(proctrialsCounter,1:size(ps2,2)) = alpha_cue; alphaTime(proctrialsCounter,size(ps2,2)+1:size(ps1,2)+7) = alpha_nav;
                         gammaTime(proctrialsCounter,1:size(ps2,2)) = gamma_cue; gammaTime(proctrialsCounter,size(ps2,2)+1:size(ps1,2)+7) = gamma_nav;
                     end

                     for r = 1:size(t1,2) % get average x- ,y- position in each time bin
                        if r == 1
                            ind = find(trialdata(:,2) < t1(r));
                            loc(r,1) = mean(trialdata(ind,3)); loc(r,2) = mean(trialdata(ind,4));  
                        elseif r > 1
                            ind = find(trialdata(:,2) > t1(r-1) & trialdata(:,2) < t1(r));
                            loc(r,1) = mean(trialdata(ind,3)); loc(r,2) = mean(trialdata(ind,4));  
                        end

                        xCoord = find(horGridBound < loc(r,1)); xCoord = xCoord(end); 
                        yCoord = find(vertGridBound > loc(r,2)); yCoord = yCoord(end);
                        loc(r,3) = xCoord + (yCoord-1)*res; % get location in each time bin

                        % Frequency bands of interest: theta (4-8 Hz), alpha (8–13 Hz), and gamma (30-120 Hz)
                        binCounter = binCounter + 1;
                        theta(binCounter,1) = mean(ps1(3:5,r)); theta(binCounter,2) = loc(r,3);
                        alpha(binCounter,1) = mean(ps1(6:8,r)); alpha(binCounter,2) = loc(r,3);
                        gamma(binCounter,1) = mean(ps1(17:62,r)); gamma(binCounter,2) = loc(r,3);
                     end

                    q = q + 1;
                    clear loc;
                catch % for trials that are too short
%                     disp(q);
                    q = q + 1;
                end
            else
                q = q + 1;
            end
        end
        binCounter = 0; q = 1; % reset counters

        % Get average power in each location across session
        for s = 1:numGrids
            thetaPower(p,s) = mean(theta(theta(:,2)==s,1)); 
            alphaPower(p,s) = mean(alpha(alpha(:,2)==s,1)); 
            gammaPower(p,s) = mean(gamma(gamma(:,2)==s,1));
        end
    end

    % Get 95th percentile (test significance using shuffling method)
    theta_locsig(1,1:numGrids) = thetaPower(1,:);
    theta_locsig(2,1:numGrids) = prctile(thetaPower(2:end,1:numGrids),2.5); 
    theta_locsig(3,1:numGrids) = prctile(thetaPower(2:end,1:numGrids),97.5);

    alpha_locsig(1,1:numGrids) = alphaPower(1,:);
    alpha_locsig(2,1:numGrids) = prctile(alphaPower(2:end,1:numGrids),2.5); 
    alpha_locsig(3,1:numGrids) = prctile(alphaPower(2:end,1:numGrids),97.5);

    gamma_locsig(1,1:numGrids) = gammaPower(1,:);
    gamma_locsig(2,1:numGrids) = prctile(gammaPower(2:end,1:numGrids),2.5); 
    gamma_locsig(3,1:numGrids) = prctile(gammaPower(2:end,1:numGrids),97.5);

    % Calculate selectivity index (SI)
    maxTheta = max(thetaPower(1,:)); sumTheta = 0;
    maxAlpha = max(alphaPower(1,:)); sumAlpha = 0;
    maxGamma = max(gammaPower(1,:)); sumGamma = 0;    
    for t = 1:25
        if isnan(thetaPower(1,t)) || isnan(alphaPower(1,t)) || isnan(gammaPower(1,t)),
            continue
        end
        sumTheta = sumTheta + thetaPower(1,t)/maxTheta;
        sumAlpha = sumAlpha + alphaPower(1,t)/maxAlpha;
        sumGamma = sumGamma + gammaPower(1,t)/maxGamma;
    end
    theta_locSI = (25-sumTheta)/24;
    alpha_locSI = (25-sumAlpha)/24;
    gamma_locSI = (25-sumGamma)/24;
    %% Poster Position
    posterPosition = [1.5,2; 4,1.5; 4,2.5; 2,3.5; 2,4.5; 4.5,4];
    posterPosition = posterPosition*res/5.0-0.5*ones(size(posterPosition))*(res/5-1);
    
    %% Plot band power over time and heat maps  
    figure('name','band power maps');

    % Theta power over time
    numWindows = size(thetaTime,2);
    thetaTime(end+1,:) = sum(thetaTime(:,1:numWindows)) ./ sum(thetaTime(:,1:numWindows) ~= 0);
    thetaTime = thetaTime(end,:);
    % Theta power heat map
    thetaPowerSpec = reshape(thetaPower(1,:),[res,res])'; % reshape spikeFreq matrix to map onto physical maze locations
    thetaPValue = anova1(thetaPowerSpec,[],'off');
%     for i = 1:size(posterPosition,1) % poster positions
%         plot(posterPosition(i,1),posterPosition(i,2),'r.','MarkerSize',40);
%     end


    % Alpha power over time
    alphaTime(end+1,:) = sum(alphaTime(:,1:numWindows)) ./ sum(alphaTime(:,1:numWindows) ~= 0);
    alphaTime = alphaTime(end,:);
    % Alpha power heat map
    alphaPowerSpec = reshape(alphaPower(1,:),[res,res])'; % reshape spikeFreq matrix to map onto physical maze locations
    alphaPValue = anova1(alphaPowerSpec,[],'off');
    
    % Gamma power over time
    gammaTime(end+1,:) = sum(gammaTime(:,1:numWindows)) ./ sum(gammaTime(:,1:numWindows) ~= 0);
    gammaTime = gammaTime(end,:);
    % Gamma power heat map
    gammaPowerSpec = reshape(gammaPower(1,:),[res,res])'; % reshape spikeFreq matrix to map onto physical maze locations
    gammaPValue = anova1(gammaPowerSpec,[],'off');
    
        
    [fp,channel,ext] = fileparts(pwd);
    [fp,array,ext] = fileparts(fp);
    [fp,session,ext] = fileparts(fp);
    [fp,day,ext] = fileparts(fp);
     
    name = 'bandfield';
    save([name '.mat'],'thetaTime','thetaPowerSpec','thetaPValue','theta_locSI','alphaTime','alphaPowerSpec','alphaPValue','alpha_locSI','gammaTime','gammaPowerSpec','gammaPValue','gamma_locSI','posterPosition');
    disp('Saved MAT')
    close all;

    clearvars -except directory rplTriggers lfpWaveform unityTriggers unityData thetaPower theta_locsig theta_locSI thetaTime alphaPower alpha_locsig alpha_locSI alphaTime gammaPower gamma_locsig gamma_locSI gammaTime

end