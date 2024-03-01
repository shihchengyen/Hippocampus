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
	for file in 1binData.csv rplparallel.mat umaze.mat
	do
		if [ ! -f $dirpath/$file ]; then
			scp hippocampus@cortex.nus.edu.sg:$line/$file $dirpath
		fi
	done
	curr=$(pwd)
	cd $dirpath
	echo $dirpath > batch.txt
        pvjob=$(qsub $curr/pvsubmit.pbs)
        qsub -W depend=afterok:$pvjob -v rad=$rad $curr/pvtrf.pbs
	echo Changing back to working dir: $curr
	cd $curr	
done
