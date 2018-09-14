disp('Submitting channel 1 to 32...')
ProcessLevel(rplsplit,'Levels','Day','SaveLevels',3,'SkipLFP','UseHPC','Channels',1:32,'HPCCmd','source ~/.bash_profile; source /etc/profile.d/rec_modules.sh; module load pbs; qsub $GITHUB_MATLAB/Hippocampus/Compiler/rplsplit/rsHPC_submit_file.txt');

pause(20)

disp(' ')
disp('Submitting channel 33 to 64...')
ProcessLevel(rplsplit,'Levels','Day','SaveLevels',3,'SkipLFP','UseHPC','Channels',33:64,'HPCCmd','source ~/.bash_profile; source /etc/profile.d/rec_modules.sh; module load pbs; qsub $GITHUB_MATLAB/Hippocampus/Compiler/rplsplit/rsHPC_submit_file.txt');

pause(20)
disp(' ')

disp('Submitting channel 65 to 96...')
ProcessLevel(rplsplit,'Levels','Day','SaveLevels',3,'SkipLFP','UseHPC','Channels',65:96,'HPCCmd','source ~/.bash_profile; source /etc/profile.d/rec_modules.sh; module load pbs; qsub $GITHUB_MATLAB/Hippocampus/Compiler/rplsplit/rsHPC_submit_file.txt');

pause(20)
disp(' ')

disp('Submitting channel 97 to 124...')
ProcessLevel(rplsplit,'Levels','Day','SaveLevels',3,'SkipLFP','UseHPC','Channels',97:124,'HPCCmd','source ~/.bash_profile; source /etc/profile.d/rec_modules.sh; module load pbs; qsub $GITHUB_MATLAB/Hippocampus/Compiler/rplsplit/rsHPC_submit_file.txt');

disp(' ')

disp('done...')

exit
