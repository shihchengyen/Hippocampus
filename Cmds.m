%% Processing data after each day of recordings
% This will split the raw Ripple files in each session into 
% rplraw and rplparallel objects. We will skip the LFP files
% as we will perform our own low-pass filtering on the raw data
ProcessLevel(rplsplit,'Levels','Day','SaveLevels',3,'SkipLFP','HPCCmd','source ~/.bash_profile; condor_submit ~/cbin/rs_submit_file.txt')
% If the HPC is out of space, we may want to split the Ripple file, and
% create the highpass and lfp files, but skip the spike sorting
ProcessLevel(rplsplit,'Levels','Day','SaveLevels',3,'SkipLFP','SkipSort','HPCCmd','source ~/.bash_profile; condor_submit ~/cbin/rs_submit_file.txt')
% split the files locally before sending them back to the Synology
ProcessLevel(rplsplit,'Levels','Day','SaveLevels',3,'SkipLFP','SkipSort')
% If we have to split some of the channels because rplsplit died halfway
% through, we can use the Channels argument
rs = rplsplit('auto','SaveLevels',3,'SkipParallel','SkipLFP','SkipSort','SkipAnalog','Channels',27)
% For longer sessions, splitting the file might take longer 
% than the 24 hour limit on the HPC, so we can
% split by arrays instead
ProcessLevel(rplsplit,'Levels','Day','SaveLevels',3,'SkipLFP','Channels',1:32,'HPCCmd','source ~/.bash_profile; condor_submit ~/cbin/rs_submit_file.txt')
% need to pause to make sure job starts before 
ProcessLevel(rplsplit,'Levels','Day','SaveLevels',3,'SkipLFP','SkipParallel','SkipAnalog','Channels',33:64,'HPCCmd','source ~/.bash_profile; condor_submit ~/cbin/rs_submit_file.txt')
ProcessLevel(rplsplit,'Levels','Day','SaveLevels',3,'SkipLFP','SkipParallel','SkipAnalog','Channels',65:96,'HPCCmd','source ~/.bash_profile; condor_submit ~/cbin/rs_submit_file.txt')
ProcessLevel(rplsplit,'Levels','Day','SaveLevels',3,'SkipLFP','SkipParallel','SkipAnalog','Channels',97:124,'HPCCmd','source ~/.bash_profile; condor_submit ~/cbin/rs_submit_file.txt')
% This will create the highpass and lfp files, but not run the spike
% sorting. Useful when debugging noise
ProcessLevel(rplsplit,'Levels','Day','Exclude',{'session01','session02','session03','session04','sessioneye'},'SaveLevels',3,'SkipLFP','SkipSort','HPCCmd','source ~/.bash_profile; condor_submit ~/cbin/rs_submit_file.txt')
% This will generate frequency plots for highpass and lfp data
ProcessLevel(nptdata,'Levels','Day','nptLevelCmd',{'Session','checkRecording'});
% Generate unitymaze objects
um = ProcessLevel(unitymaze,'Levels','Days','redo','Include',{'20180816','20180817','2018082','201809'});
InspectGUI(um)
figure
for i = 1:9
	plot(um,i)
    saveas(gcf,['um' num2str(i) '.png'])
end

ProcessLevel(rplsplit,'Levels','Day','SaveLevels',3,'Test','SkipLFP','SkipRaw','SkipAnalog','SkipParallel','HPCCmd','source ~/.bash_profile; condor_submit ~/cbin/rs_submit_file.txt')
ProcessLevel(rplsplit,'Levels','Session','SaveLevels',3,'SkipSort','SkipLFP','SkipParallel','Channels',91:124,'HPCCmd','source ~/.bash_profile; condor_submit ~/cbin/rs_submit_file.txt')

ProcessLevel(nptdata,'Levels','Day','Exclude',{'analog'},'nptLevelCmd',{'Channel','system(''source ~/.bash_profile; condor_submit ~/cbin/hplfp_submit_file.txt'')'})
ProcessLevel(nptdata,'Levels','Session','Exclude',{'analog'},'nptLevelCmd',{'Channel','sortchannel'})

ProcessLevel(nptdata,'Levels','Session','nptLevelCmd',{'Session','submitSort(''HPC'')'})

% generate rpllfp objects in eyesession directories
ProcessLevel(nptdata,'Levels','Session','nptLevelCmd',{'Channel','rpllfp(''auto'',''save'');'})

% Generate plots by arrays of the LFPs
% Break long trials into plots of 3.5 s
ProcessLevel(nptdata,'Levels','Session','nptLevelCmd',{'Array', ...
	'cwd = pwd; cd ../analog; vma = ProcessLevel(vmlfp,''Levels'',''Array'',''SaveLevels'',2); cd(cwd); vmb = ProcessLevel(vmlfp,''Levels'',''Array''); vmc = vmb + vma; sindex(1).type = ''.''; sindex(1).subs = ''data''; sindex(2).type = ''.''; sindex(2).subs = ''timeStamps''; tStamps = subsref(vmc,sindex); ntrials = size(tStamps,1); a = nptSubplot(1,1); h = gcf; h.PaperOrientation = ''landscape''; h.PaperType = ''a4''; h.Visible = ''off''; for i = 1:ntrials; plot(vmc,i); linkaxes(get(gcf,''Children''),''x''); title(pwd); xmaxi = tStamps(i,end)-tStamps(i,1)+0.5; ntplots = ceil(xmaxi/3.5); xlimits = -500:3500:((ntplots*3500)-500); for j = 1:ntplots; xlim([xlimits(j) xlimits(j+1)]); print(''-dpdf'',''-fillpage'',''-noui'',[getDataOrder(''ShortName'') sprintf(''LFP_Analog_t%03d'',i) sprintf(''-%02d.pdf'',j)]); end; subplot(a(1),a(2),1); end'})

ProcessLevel(nptdata,'Levels','Session','nptLevelCmd',{'Array', ...
	'ng = nptgroup(''auto'',''Levels'',''Array'',''LevelObject'',''Channel''); cwd = pwd; dirlist = dir(''channel*''); cd(dirlist(1).name); vp = vmlfp(''auto''); ntrials = vp.number; sindex(1).type = ''.''; sindex(1).subs = ''data''; sindex(2).type = ''.''; sindex(2).subs = ''timeStamps''; tStamps = subsref(vp,sindex); cd(cwd); a = nptSubplot(1,1); h = gcf; h.PaperOrientation = ''landscape''; h.PaperType = ''a4''; for i = 1:ntrials; plot(ng,i,''Object'',{''vmlfp''},''GroupPlotSep'',''Vertical''); linkaxes(get(gcf,''Children''),''x''); title(pwd); xmaxi = tStamps(i,end)-tStamps(i,1)+0.5; ntplots = ceil(xmaxi/3.5); xlimits = -500:3500:((ntplots*3500)-500); for j = 1:ntplots; xlim([xlimits(j) xlimits(j+1)]); print(''-dpdf'',''-fillpage'',''-noui'',[getDataOrder(''ShortName'') sprintf(''LFPt%03d'',i) sprintf(''-%02d.pdf'',j)]); end; subplot(a(1),a(2),1); end'})

% Generate plots by arrays of the highpass data
% Break long trials into plots of 3.5 s
ProcessLevel(nptdata,'Levels','Session','nptLevelCmd',{'Array', ...
	'ng = nptgroup(''auto'',''Levels'',''Array'',''LevelObject'',''Channel''); cwd = pwd; dirlist = dir(''channel*''); cd(dirlist(1).name); vp = vmhighpass(''auto''); ntrials = vp.number; sindex(1).type = ''.''; sindex(1).subs = ''data''; sindex(2).type = ''.''; sindex(2).subs = ''timeStamps''; tStamps = subsref(vp,sindex); cd(cwd); a = nptSubplot(1,1); h = gcf; h.PaperOrientation = ''landscape''; h.PaperType = ''a4''; for i = 1:ntrials; plot(ng,i,''Object'',{''vmhighpass''},''GroupPlotSep'',''Vertical''); linkaxes(get(gcf,''Children''),''x''); title(pwd); xmaxi = tStamps(i,end)-tStamps(i,1)+0.5; ntplots = ceil(xmaxi/3.5); xlimits = -500:3500:((ntplots*3500)-500); for j = 1:ntplots; xlim([xlimits(j) xlimits(j+1)]); print(''-dpdf'',''-fillpage'',''-noui'',[getDataOrder(''ShortName'') sprintf(''HPt%03d'',i) sprintf(''-%02d.pdf'',j)]); end; subplot(a(1),a(2),1); end'})

ProcessLevel(rplsplit,'Levels','Session','SaveLevels',3,'SkipLFP')

% check 50 Hz and high-frequency noise
checkRecording

% check spike sorting
vs = ProcessLevel(viewsort,'Levels','Session');

% You can plot the waveforms by channel:
InspectGUI(vs)

% You can also plot the waveforms by array:
InspectGUI(vs,'Array')

% To plot the waveforms by array using the GMR geometry, you can do:
hd = ProcessLevel(hidata,'Levels','Session','FileName','hmmsort.mat');
InspectGUI(hd,'Array','UseGMR','Objects',{'viewsort',{},{}})

% You can plot waveforms, see highpass data and save spiketrains by channel:
InspectGUI(vs,'Cmds','InspectGUI(vmhighpass(''auto''),''LoadSort''); pause; system(''source ~/.bash_profile; /opt/data2/anaconda2/bin/python ~/Dropbox/Work/Matlab/hmmsort/hmmsort/create_spiketrains.py'')')
% Add vmplacecell plot into the above
InspectGUI(vs,'Cmds','InspectGUI(vmhighpass(''auto''),''LoadSort''); vpc=vmplacecell(''auto'',''GridSteps'',10,''SaveLevels'',2,''MinTrials'',0);InspectGUI(vpc,''addObjs'',{vpc},''optArgs'',{{},{''Errorbar''}},''SP'',[2 1]); pause; system(''source ~/.bash_profile; /opt/data2/anaconda2/bin/python ~/Dropbox/Work/Matlab/hmmsort/hmmsort/create_spiketrains.py'')')

% To create the viewsort object using saved spiketrains, then plot as before:
vss = ProcessLevel(viewsort, 'Levels', 'Session', 'Saved')
InspectGUI(vss,'Array')

% This will use rplraw to first create rplhighpass, then use rplparallel
% to create vmhighpass
ProcessLevel(vmhighpass,'Levels','Session','SaveLevels',2)
% This will use rplraw to first create rpllfp, then use rplparallel
% to create vmlfp
ProcessLevel(vmlfp,'Levels','Session','SaveLevels',2)
ProcessLevel(nptdata,'Levels','Session','nptLevelCmd',{'Array', ...
	'ng = nptgroup(''auto'',''Levels'',''Array'',''LevelObject'',''Channel''); cwd = pwd; dirlist = dir(''channel*''); cd(dirlist(1).name); vp = vmlfp(''auto''); ntrials = vp.number; cd(cwd); a = nptSubplot(1,1); h = gcf; h.PaperOrientation = ''landscape''; h.PaperType = ''a4''; for i = 1:ntrials; plot(ng,i,''Object'',{''vmlfp''},''GroupPlotSep'',''Vertical''); title(pwd); print(''-dpdf'',''-fillpage'',''-noui'',[getDataOrder(''ShortName'') sprintf(''LFPt%03d.pdf'',i)]); subplot(a(1),a(2),1); end'})

ProcessLevel(vmhighpass,'Levels','Session','SaveLevels',2)
ProcessLevel(nptdata,'Levels','Session','Exclude',{'array01'},'nptLevelCmd',{'Array', ...
	'ng = nptgroup(''auto'',''Levels'',''Array'',''LevelObject'',''Channel''); cwd = pwd; dirlist = dir(''channel*''); cd(dirlist(1).name); vp = vmhighpass(''auto''); ntrials = vp.number; cd(cwd); a = nptSubplot(1,1); h = gcf; h.PaperOrientation = ''landscape''; h.PaperType = ''a4''; for i = 1:ntrials; plot(ng,i,''Object'',{''vmhighpass''},''GroupPlotSep'',''Vertical''); title(pwd); print(''-dpdf'',''-fillpage'',''-noui'',[getDataOrder(''ShortName'') sprintf(''HPt%03d.pdf'',i)]); subplot(a(1),a(2),1); end'})

ProcessLevel(nptdata,'Levels','Array','nptLevelCmd',{'Array', ...
	'ng = nptgroup(''auto'',''Levels'',''Array'',''LevelObject'',''Channel''); cwd = pwd; dirlist = dir(''channel*''); cd(dirlist(1).name); vp = vmhighpass(''auto''); ntrials = vp.number; cd(cwd); a = nptSubplot(1,1); h = gcf; h.PaperOrientation = ''landscape''; h.PaperType = ''a4''; for i = 1:ntrials; plot(ng,i,''Object'',{''vmhighpass''},''GroupPlotSep'',''Vertical''); title(pwd); print(''-dpdf'',''-fillpage'',''-noui'',[getDataOrder(''ShortName'') sprintf(''HPt%03d.pdf'',i)]); subplot(a(1),a(2),1); end'})

ng = nptgroup('auto','Levels','Array','LevelObject','Channel');
ntrials = 82;
a = nptSubplot(1,1); h = gcf; h.PaperOrientation = 'landscape'; h.PaperType = 'a4';
for i = 2:ntrials
    plot(ng,i,'Object',{'vmhighpass'},'GroupPlotSep','Vertical'); 
    title(pwd); 
    print('-dpdf','-fillpage','-noui',[getDataOrder('ShortName') sprintf('HPt%03d.pdf',i)]); 
    subplot(a(1),a(2),1); 
end

% check spike sorting
h = rplhighpass('auto');
load g095c01_spiketrain.mat
panGUI
plot(rh.data.analogTime(1:100000)*1000,rh.data.analogData(1:100000),'.-')
hold on
stem(rh.data.analogTime(timestamps(1:100))*1000,repmat(100,1,100))

% Create nptgroup object containing all channels in a session
cd('session01')
cd('array01')
ng1 = nptgroup('auto','Levels','Array','LevelObject','Channel')
cd('../array02')
ng2 = nptgroup('auto','Levels','Array','LevelObject','Channel')
cd('../array03')
ng3 = nptgroup('auto','Levels','Array','LevelObject','Channel')
cd('../array04')
ng4 = nptgroup('auto','Levels','Array','LevelObject','Channel')
InspectGUI(ng1,'Object',{'vmlfp'},'GroupPlotSep','Vertical')
InspectGUI(ng2,'Object',{'vmlfp'},'GroupPlotSep','Vertical')
InspectGUI(ng3,'Object',{'vmlfp'},'GroupPlotSep','Vertical')
InspectGUI(ng4,'Object',{'vmlfp'},'GroupPlotSep','Vertical')

% Inspect all vmhighpass objects session by session
InspectGUI(ng27_1,'Object',{'vmhighpass'},'GroupPlotSep','Vertical')
% Inspect all vmlfp objects
InspectGUI(ng27_1,'Object',{'vmlfp'},'GroupPlotSep','Vertical')
InspectGUI(ng27_2,'Object',{'vmlfp'},'GroupPlotSep','Vertical')
InspectGUI(ng27_3,'Object',{'vmlfp'},'GroupPlotSep','Vertical')
InspectGUI(ng25_1,'Object',{'vmlfp',{'OldMarkerFormat2'}},'GroupPlotSep','Vertical')
InspectGUI(ng25_2,'Object',{'vmlfp',{'OldMarkerFormat2'}},'GroupPlotSep','Vertical')
InspectGUI(ng25_3,'Object',{'vmlfp',{'OldMarkerFormat2'}},'GroupPlotSep','Vertical')
InspectGUI(ng25_4,'Object',{'vmlfp',{'OldMarkerFormat2'}},'GroupPlotSep','Vertical')
InspectGUI(ng25_5,'Object',{'vmlfp',{'OldMarkerFormat2'}},'GroupPlotSep','Vertical')
InspectGUI(ng18_1_1,'Object',{'vmlfp',{'OldMarkerFormat'},{'LowpassFreqs',[0.3 250]}},'GroupPlotSep','Vertical')


%% Check Ripple file using Neuroshare functions
[nsStatus,hFile] = ns_OpenFile('171219_Block1.ns5')
entityNum = 256;
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
rplsplit('auto','SaveLevels',3,'SkipLFP');
rplsplit('auto','SaveLevels',2,'SkipRaw','LowpassFreqs',[3 250]);
rplsplit('auto','redo','SaveLevels',2,'SkipLFP','SkipParallel','Channels',114);
ProcessLevel(rplsplit,'Levels','Day','SaveLevels',2,'LevelObject','Session','redo')
ProcessLevel(rplsplit,'Levels','Days','SaveLevels',3,'LevelObject','Session','SkipRaw','SkipLFP','SkipParallel')


%% create LFP from raw data
% run from channel data containing rplraw.mat
rp = rpllfp('auto','save')
vp = vmlfp('auto','LowpassFreqs',[3 250])

% create vmlfp objects in each channel
ProcessLevel(vmlfp,'Levels','session','save')
ProcessLevel(vmhighpass,'Levels','session','SaveLevels',2)
ProcessLevel(vmhighpass,'Levels','Day','SaveLevels',2)

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
InspectGUI(ng,'Object',{'vmlfp',{'PreTrial',750},{'LowpassFreqs',[3 250]}},'GroupPlotSep','Vertical')
InspectGUI(ng,'Object',{'vmlfp',{'FreqPlot','PlotAllData','FreqLims',[0 150]}}, ...
	'GroupPlotSep','Vertical')


cd('channel022')
vp = vmlfp('auto')
% Plot LFP and frequency spectrum
InspectGUI(vp,'addObjs',{vp},'optArgs',{{},{'FreqPlot','FreqLims',[0 150]}},'SP',[2 1])
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

[nptThesisCells,vpc20] = ProcessDirs(nptThesisCells,'Object','vmplacecell','save','redo','GridSteps',20);
ProcessDirs(nptThesisSessions,'Object','unitymaze','redo','save','GridSteps',10)
[nptThesisCells,vpc10] = ProcessDirs(nptThesisCells,'Object','vmplacecell','save','redo','GridSteps',10);
ProcessDirs(nptThesisSessions,'Object','unitymaze','redo','save','GridSteps',5)
[nptThesisCells,vpc5] = ProcessDirs(nptThesisCells,'Object','vmplacecell','save','redo','GridSteps',5);

