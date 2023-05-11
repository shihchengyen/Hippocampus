#!/bin/bash

rad=$1
cat ../sessions.txt | while read line
do
	temp="${line##*Data}/session01"
	echo Reading from sessions.txt: $line
	line+="/session01"
	dirpath=~/hpctmp/Data$temp
	mkdir -p $dirpath
	echo hpc session directory: $dirpath
	if [ ! -f $dirpath/unityfile.mat ]; then
		scp -P 8398 hippocampus@cortex.nus.edu.sg:$line/unityfile.mat $dirpath 
	fi
	if [ ! -f $dirpath/eyelink.mat ]; then
		scp -P 8398 hippocampus@cortex.nus.edu.sg:$line/eyelink.mat $dirpath
	fi
	curr=$(pwd)
	cd $dirpath
	echo $dirpath > batch.txt
		rcjob=$(qsub -v rad=$rad $curr/rcsubmit.pbs)
		qsub -W depend=afterok:$rcjob -v rad=$rad $curr/rctrf.pbs
	echo Changing back to working dir: $curr
	cd $curr	
done

