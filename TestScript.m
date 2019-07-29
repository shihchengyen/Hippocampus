ProcessLevel(nptdata,'Levels','Days','Include',{'20181102'},'Exclude',{'sessioneye','session_merged','sessiontest'},'nptLevelCmd',{'Channel','vmlfp(''auto'',''redo'');'});
ProcessLevel(nptdata,'Levels','Days','Include',{'20181102'},'Exclude',{'sessioneye','session_merged','sessiontest'},'nptLevelCmd',{'Session','unitymaze(''auto'',''redo'');'});
ProcessLevel(vmplacecell,'Levels','Days','Include',{'20181102'});
ProcessLevel(nptdata,'Levels','Days','Include',{'20181102'},'Exclude',{'sessioneye','session_merged','sessiontest'},'nptLevelCmd',{'Session','eyelink(''auto'',''redo'');'});
cd 20181102/session01
um = unitymaze('auto','GridSteps',40);
el = eyelink('auto');
InspectGUI(um,'addObjs',{el},'optArgs',{{'Trial'},{'Trial'}},'SP',[2 1])
vs = ProcessLevel(viewsort,'Levels','Session');
InspectGUI(vs)
InspectGUI(vs,'Array')
InspectGUI(vs,'Cmds','InspectGUI(vmhighpass(''auto''),''LoadSort''); vpc=vmplacecell(''auto'',''GridSteps'',10,''SaveLevels'',2);InspectGUI(vpc,''addObjs'',{vpc},''optArgs'',{{},{''Errorbar''}},''SP'',[2 1])')
hd = ProcessLevel(hidata,'Levels','Session','FileName','hmmsort.mat');
InspectGUI(hd,'Array','UseGMR','Objects',{'viewsort',{},{}})
cd array01
ProcessLevel(nptdata,'Levels','Array','nptLevelCmd',{'Channel','vmhighpass(''auto'',''redo'');'});
cd channel010
vl = vmlfp('auto');
InspectGUI(um,'addObjs',{el,vl},'optArgs',{{'Trial'},{'Trial'},{}},'SP',[3 1])
cd ../../../../
% unityfile
% raycast
