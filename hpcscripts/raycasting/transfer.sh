#!/bin/bash

curr=$(pwd)
rad=$1
cat sessions.txt | while read line
do
	temp=${line##*Data}
	echo Transferring files from: ~/hpctmp/Data$temp
	cd ~/hpctmp/Data$temp
	scp -P 8398 binData.hdf hippocampus@cortex.nus.edu.sg:$line/${rad}binData.hdf 
	scp -P 8398 unityfile_eyelink.csv logs.txt VirtualMazeBatchLog.txt hippocampus@cortex.nus.edu.sg:$line/ 
done
cd $curr

