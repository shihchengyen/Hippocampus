disp('Submitting channel 1 to 32...')
ProcessLevel(rplsplit,'skipCheckingMarkers','Levels','Day','SaveLevels',3,'SkipLFP','UseHPC','Channels',1:32,'skipCheckingRplsplit','HPCCmd','source ~/.bash_profile; qsub $GITHUB_MATLAB/Hippocampus/Compiler/rplsplit/rsHPC_submit_file.txt');

pause(20)

disp(' ')
disp('Submitting channel 33 to 64...')
ProcessLevel(rplsplit,'skipCheckingMarkers','Levels','Day','SaveLevels',3,'SkipLFP','UseHPC','Channels',33:64,'skipCheckingRplsplit','HPCCmd','source ~/.bash_profile; qsub $GITHUB_MATLAB/Hippocampus/Compiler/rplsplit/rsHPC_submit_file.txt');

pause(20)
disp(' ')

disp('Submitting channel 65 to 96...')
ProcessLevel(rplsplit,'skipCheckingMarkers','Levels','Day','SaveLevels',3,'SkipLFP','UseHPC','Channels',65:96,'skipCheckingRplsplit','HPCCmd','source ~/.bash_profile; qsub $GITHUB_MATLAB/Hippocampus/Compiler/rplsplit/rsHPC_submit_file.txt');

pause(20)
disp(' ')

disp('Submitting channel 97 to 124...')
ProcessLevel(rplsplit,'skipCheckingMarkers','Levels','Day','SaveLevels',3,'SkipLFP','UseHPC','Channels',97:124,'skipCheckingRplsplit','HPCCmd','source ~/.bash_profile; qsub $GITHUB_MATLAB/Hippocampus/Compiler/rplsplit/rsHPC_submit_file.txt');
disp(' ')

disp('done...')

exit
