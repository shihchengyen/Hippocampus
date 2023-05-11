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
	if [ ! -f $dirpath/rplparallel.mat ]; then
    	scp -P 8398 hippocampus@cortex.nus.edu.sg:$line/rplparallel.mat $dirpath
	fi
	if [ ! -f $dirpath/unityfile.mat ]; then
		scp -P 8398 hippocampus@cortex.nus.edu.sg:$line/unityfile.mat $dirpath 
	fi
	if [ ! -f $dirpath/eyelink.mat ]; then
		scp -P 8398 hippocampus@cortex.nus.edu.sg:$line/eyelink.mat $dirpath
	fi
	if [ ! -f $dirpath/umaze.mat ]; then
		scp -P 8398 hippocampus@cortex.nus.edu.sg:$line/umaze.mat $dirpath
	fi
	curr=$(pwd)
	cd $dirpath
	echo $dirpath > batch.txt
		rcjob=$(qsub -v rad=$rad $curr/rcsubmit.pbs)
		pvjob=$(qsub -W depend=afterok:$rcjob $curr/pvsubmit.pbs)
		qsub -W depend=afterok:$pvjob -v rad=$rad $curr/rcpvtrf.pbs
	echo Changing back to working dir: $curr
	cd $curr	
done

