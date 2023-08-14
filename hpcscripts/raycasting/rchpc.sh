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
	for file in unityfile.mat eyelink.mat
	do
		if [ ! -f $dirpath/$file ]; then
			scp hippocampus@cortex.nus.edu.sg:$line/$file $dirpath 
		fi
	done
	curr=$(pwd)
	cd $dirpath
	echo $dirpath > batch.txt
		rcjob=$(qsub -v rad=$rad $curr/rcsubmit.pbs)
		qsub -W depend=afterok:$rcjob -v rad=$rad $curr/rctrf.pbs
	echo Changing back to working dir: $curr
	cd $curr	
done

