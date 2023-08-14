#!/bin/bash

curr=$(pwd)
rad=$1
cat sessions.txt | while read line
do
	temp=${line##*Data}
	echo Transferring files from: ~/hpctmp/Data$temp
	cd ~/hpctmp/Data$temp
	scp vmpv.mat hippocampus@cortex.nus.edu.sg:$line/${rad}vmpv.mat 
done
cd $curr

