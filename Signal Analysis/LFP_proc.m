raw_LFP = analogInputData_uV;

% Remove data packet loss distortion
raw_LFP(raw_LFP < -1000) = 0;

%% Method from 'Sharp wave ripples during visual exploration in the primate hippocampus, Leonard et al. 2015
% Band-pass filter (100-250Hz)
bandPass_LFP_meth1 = ft_preproc_bandpassfilter(raw_LFP,1000,[100 250]);
% Transform to z-scores
zscore_LFP_meth1 = zscore(bandPass_LFP_meth1);
% Rectify signal
rectified_LFP_meth1 = abs(zscore_LFP_meth1);
% Low-pass filter (1-20Hz)
bandPass2_LFP_meth1 = ft_preproc_bandpassfilter(rectified_LFP_meth1,1000,[1 20]);
% Identify ripples, threshold crossings 3 SDs above the mean, with a minimum duration of 50ms, beginning and ending at 1 SD

%% Method from 'Midline thalamic neurons are differentially engaged during HC network oscillations', Lara-Vasquez et al. 2016
% Check that LFP is downsampled to 1kHz
% Band-pass filter (100-200Hz) using zero phase shift non-causal FIR filter
% with 0.5Hz roll-off
bandPass_LFP_meth2 = ft_preproc_bandpassfilter(raw_LFP,1000,[100 200],6,'fir');
% Rectify signal
rectified_LFP_meth2 = abs(bandPass_LFP_meth2);
% Low-pass filter at 20Hz with 4th order Butterworth filter
lowPass_LFP_meth2 = ft_preproc_lowpassfilter(rectified_LFP_meth2,1000,20,4,'but');
% z-score normalize using mean and SD of whole signal in the time domain
zscore_LFP_meth2 = zscore(lowPass_LFP_meth2);
% Epochs during which normalized signal exceeded 3.5SD threshold were
% considered ripple events

%% Plot figures
xMin = 100000;
xMax = 101000;

figure()
plot(raw_LFP); 
xlim([xMin xMax]); % plot 10 secs, 100th-110th second
xlabel('Time in milliseconds'); ylabel('microVolts');
title('Raw LFP'); 

figure()
subplot(2,2,1)
plot(bandPass_LFP_meth1);
xlim([xMin xMax]); % plot 10 secs, 100th-110th second
title('Band-pass filter (100-250Hz)');

subplot(2,2,2)
plot(zscore_LFP_meth1);
xlim([xMin xMax]); % plot 10 secs, 100th-110th second
title('z-transform signal');

subplot(2,2,3)
plot(rectified_LFP_meth1);
xlim([xMin xMax]); % plot 10 secs, 100th-110th second
title('Rectified signal');

subplot(2,2,4)
plot(bandPass2_LFP_meth1);
xlim([xMin xMax]); % plot 10 secs, 100th-110th second
title('Band-pass filter (1-20Hz)'); 

figure()
subplot(2,2,1)
plot(bandPass_LFP_meth2);
xlim([xMin xMax]); % plot 10 secs, 100th-110th second
title('Band-pass finite impulse filter (100-200Hz)'); 

subplot(2,2,2)
plot(rectified_LFP_meth2);
xlim([xMin xMax]); % plot 10 secs, 100th-110th second
title('Rectified signal'); 

subplot(2,2,3)
plot(lowPass_LFP_meth2);
xlim([xMin xMax]); % plot 10 secs, 100th-110th second
title('Low-pass 4th order Butterworth filter (20Hz)'); 

subplot(2,2,4)
plot(zscore_LFP_meth2);
xlim([xMin xMax]); % plot 10 secs, 100th-110th second
title('z-transform signal'); 

% time in x, frequency in y (time-frequency analysis), spectrogram, evoked
% minus spontaneous spectrum to take care of baseline (but this may be more
% useful for comparing conditions rather than baseline because ripple
% frequencies may be present at baseline too)