% run in session directory to create rplraw and rpllfp objects in each channel
% if we want to just create rpllfp and not rplraw objects
rplsplit('auto','SaveLevels',2,'SkipRaw');

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
