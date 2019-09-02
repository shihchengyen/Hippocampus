#!/bin/sh

cwd=`pwd`; 
for i in `find . -name "channel*"`; do 
	echo $i; 
	cd $i; 
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
	cd $cwd; 
done
