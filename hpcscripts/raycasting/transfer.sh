#!/bin/bash

curr=$(pwd)
rad=$1
cat sessions.txt | while read line
do
	temp=${line##*picasso-misc}
	echo Transferring files from: ~/hpctmp/Data/picasso-misc$temp
	cd ~/hpctmp/Data/picasso-misc$temp
	scp ${rad}binData.csv hippocampus@cortex.nus.edu.sg:/picasso-misc-folder-link$temp/${rad}binData.csv
	scp unityfile_eyelink.csv logs.txt VirtualMazeBatchLog.txt hippocampus@cortex.nus.edu.sg:/picasso-misc-folder-link$temp/ 
done
cd $curr

