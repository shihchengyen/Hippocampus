#!/bin/sh

#dataCheckList.txt is saved under ~/Documents
listPath=~/Documents/dataCheckList.txt

while read line; do
	year="/Volumes"
	index="${line%%$year*}"
	sessionDir=${line:${#index}}

      	if ! [[ $line == *"load"* ]]; then
		cd $sessionDir
		if [ -f rplraw.mat ] && [ -f rpllfp.mat ] && [ -f rplhighpass.mat ]
		then
               		condor_submit $GITHUB_MATLAB/Hippocampus/Compiler/hplfp/resaveAll_submit_file.txt; 
		elif [ -f rplraw.mat ] && [ -f rpllfp.mat ]
   		then
			condor_submit $GITHUB_MATLAB/Hippocampus/Compiler/hplfp/resaveRawAndLfp_submit_file.txt;
		elif [ -f rplraw.mat ]
		then
			condor_submit $GITHUB_MATLAB/Hippocampus/Compiler/hplfp/resaveRaw_submit_file.txt;
		fi	
	fi
done < $listPath

