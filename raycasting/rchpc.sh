#!/bin/bash

cat sessions.txt | while read line
do
	temp=${line##*Data}
	echo $line
	echo ~/hpctmp/Data$temp
	scp -P 8398 hippocampus@cortex.nus.edu.sg:$line/unityfile.mat ~/hpctmp/Data$temp 
	scp -P 8398 hippocampus@cortex.nus.edu.sg:$line/eyelink.mat ~/hpctmp/Data$temp
	curr=$(pwd)
	echo $curr
	cd ~/hpctmp/Data$temp
	echo ~/hpctmp/Data$temp > batch.txt
       	qsub $curr/rcsubmit.pbs	
	cd $curr	
done
