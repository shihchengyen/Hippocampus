%% Check Ripple file using Neuroshare functions
[nsStatus,hFile] = ns_OpenFile('170918.ns5');
entityNum = 279;
[nsStatus,entityInfo] = ns_GetEntityInfo(hFile,entityNum)
[nsStatus,analogInfo] = ns_GetAnalogInfo(hFile,entityNum)
[nsStatus,~,analogData] = ns_GetAnalogData(hFile,entityNum,1,entityInfo.ItemCount);
ns_status = ns_CloseFile(hFile);

%% Inspect Ripple file
% run from session directory containing ns5 file
rv = rplview('auto')
InspectGUI(rv)

%% run in session directory to create rplraw and rpllfp objects in each channel
% if we want to just create rpllfp and not rplraw objects
rplsplit('auto','SaveLevels',2,'SkipRaw');
rplsplit('auto','SaveLevels',2,'SkipLFP','SkipParallel','Channels',1:10);

%% create LFP from raw data
% run from channel data containing rplraw.mat
rp = rpllfp('auto','save')
vp = vmlfp('auto')

% create vmlfp objects in each channel
ProcessLevel(vmlfp,'Levels','session','save')

% change directory to one of the arrays, e.g. array01
cd('array01')
% create objects
cd('channel001')
vp1 = vmlfp('auto');

% look at Channel 1 trial-by-trial
InspectGUI(vp1)
% increase the pre-trial window to 1000 ms
InspectGUI(vp1,'PreTrial',1000)
% plot all trials as a heat map
plot(vp1)
% normalize each trial before plotting heatmap
plot(vp1,'NormalizeTrial')

cd('../channel002');
vp2 = vmlfp('auto');

% look at Channels 1 and 2 together, with the subplots arranged as [2 1]
InspectGUI(vp1,'addObjs',{vp2},'SP',[2 1])

% look at all channels on array 1
% change to array01 directory
cd ..
% create nptgroup object
ng = nptgroup('auto','CellName','channel0*');
% plot all channels in array 1 (arranged vertically) using a pre-trial window of 750 ms
InspectGUI(ng,'Object',{'vmlfp',{'PreTrial',750}},'GroupPlotSep','Vertical')
InspectGUI(ng,'Object',{'vmlfp',{'FreqPlot','PlotAllData','FreqLims',[0 150]}}, ...
	'GroupPlotSep','Vertical')


cd('channel022')
vp = vmlfp('auto')
% Plot LFP and frequency spectrum
InspectGUI(vp,'addObjs',{vp},'optArgs',{{},{'FreqPlot'}},'SP',[2 1])
% Plot LFP, cleaned LFP, frequency spectrum, and frequency spectrum for 
% cleaned LFP
InspectGUI(vp,'addObjs',{vp,vp,vp},'optArgs',{{},{'RemoveLineNoise',50}, ...
	{'FreqPlot'},{'FreqPlot','RemoveLineNoise',50}},'SP',[4 1])
% Limit the frequency axis to 0 to 150
InspectGUI(vp,'addObjs',{vp,vp,vp},'optArgs',{{},{'RemoveLineNoise',50}, ...
	{'FreqPlot','FreqLims',[0 150]}, {'FreqPlot','RemoveLineNoise',50, ...
	'FreqLims',[0 150]}},'SP',[4 1])
% plot cleaned LFP and time-frequency plot of cleaned LFP
InspectGUI(vp,'addObjs',{vp},'optArgs',{{'RemoveLineNoise',50}, ...
	{'TFfft','RemoveLineNoise',50}},'SP',[2 1])
% plot averaged time-frequency results
InspectGUI(vp,'TFfft','RemoveLineNoise',50,'PlotAllData')
