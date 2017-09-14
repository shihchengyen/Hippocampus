% run in session directory to create rplraw and rpllfp objects in each channel
rplsplit('auto','SaveLevels',2);

% if we want to just create rpllfp and not rplraw objects
rplsplit('auto','redo','SaveLevels',2,'SkipRaw');

% create vmlfp objects in each channel
ProcessLevel(vmlfp,'Levels','session','save')

% create objects
cd('channel001')
vp1 = vmlfp('auto');
cd('../channel002');
vp2 = vmlfp('auto');
cd('../channel003');
vp3 = vmlfp('auto');
cd('../channel004');
vp4 = vmlfp('auto');
cd('../channel005');
vp5 = vmlfp('auto');
cd('../channel006');
vp6 = vmlfp('auto');
cd('../channel007');
vp7 = vmlfp('auto');
cd('../channel008');
vp8 = vmlfp('auto');
cd('../channel009');
vp9 = vmlfp('auto');
cd('../channel010');
vp10 = vmlfp('auto');
cd('../channel011');
vp11 = vmlfp('auto');
cd('../channel012');
vp12 = vmlfp('auto');
cd('../channel013');
vp13 = vmlfp('auto');
cd('../channel014');
vp14 = vmlfp('auto');
cd('../channel015');
vp15 = vmlfp('auto');
cd('../channel016');
vp16 = vmlfp('auto');
cd('../channel017');
vp17 = vmlfp('auto');
cd('../channel018');
vp18 = vmlfp('auto');
cd('../channel019');
vp19 = vmlfp('auto');
cd('../channel020');
vp20 = vmlfp('auto');
cd('../channel021');
vp21 = vmlfp('auto');
cd('../channel022');
vp22 = vmlfp('auto');
cd('../channel023');
vp23 = vmlfp('auto');
cd('../channel024');
vp24 = vmlfp('auto');
cd('../channel025');
vp25 = vmlfp('auto');
cd('../channel026');
vp26 = vmlfp('auto');
cd('../channel027');
vp27 = vmlfp('auto');
cd('../channel028');
vp28 = vmlfp('auto');
cd('../channel029');
vp29 = vmlfp('auto');
cd('../channel030');
vp30 = vmlfp('auto');
cd('../channel031');
vp31 = vmlfp('auto');
cd('../channel032');
vp32 = vmlfp('auto');

InspectGUI(vp1,'addObjs',{vp2,vp3,vp4,vp5,vp6,vp7,vp8,vp9,vp10,vp11,vp12, ...
	vp13,vp14,vp15,vp16,vp17,vp18,vp19,vp20,vp21,vp22,vp23,vp24,vp25,vp26,...
	vp27,vp28,vp29,vp30,vp31,vp32})
	