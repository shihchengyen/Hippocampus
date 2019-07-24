ProcessLevel(nptdata,'Levels','Days','Include',{'20181102'},'Exclude',{'sessioneye'},'nptLevelCmd',{'Channel','vmlfp(''auto'',''redo'')'});
ProcessLevel(nptdata,'Levels','Days','Include',{'20181102'},'Exclude',{'sessioneye'},'nptLevelCmd',{'Channel','vmhighpass(''auto'',''redo'')'});
ProcessLevel(nptdata,'Levels','Days','Include',{'20181102'},'Exclude',{'sessioneye'},'nptLevelCmd',{'Session','unitymaze(''auto'',''redo'')'});
vc = ProcessLevel(vmplacecell,'Levels','Days','Include',{'20181102'});
ProcessLevel(nptdata,'Levels','Days','Include',{'20181102'},'nptLevelCmd',{'Session','eyelink(''auto'',''redo'')'});
cd 20181102/session01
um = unitymaze('auto','GridSteps',40);
el = eyelink('auto');
InsepctGUI(um,'addObjs',{el},'optArgs',{{'Trial'},{'Trial'}},'SP',[2 1])
cd ../..
% unityfile
% raycast
