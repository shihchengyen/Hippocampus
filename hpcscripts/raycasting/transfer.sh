#!/bin/bash

curr=$(pwd)
rad=$1
cat sessions.txt | while read line
do
	temp=${line##*Data}
	echo Transferring files from: ~/hpctmp/Data$temp

	cd ~/hpctmp/Data$temp
	scp ${rad}binData.csv hippocampus@cortex.nus.edu.sg:$line/${rad}binData.csv
	scp unityfile_eyelink.csv logs.txt VirtualMazeBatchLog.txt hippocampus@cortex.nus.edu.sg:$line/ 
done
cd $curr

