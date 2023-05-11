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
    if [ ! -f $dirpath/binData.hdf ]; then
        scp -P 8398 hippocampus@cortex.nus.edu.sg:$line/binData.hdf $dirpath
    fi
	if [ ! -f $dirpath/rplparallel.mat ]; then
    	scp -P 8398 hippocampus@cortex.nus.edu.sg:$line/rplparallel.mat $dirpath
	fi
	if [ ! -f $dirpath/umaze.mat ]; then
		scp -P 8398 hippocampus@cortex.nus.edu.sg:$line/umaze.mat $dirpath
	fi
	curr=$(pwd)
	cd $dirpath
	echo $dirpath > batch.txt
        pvjob=$(qsub $curr/pvsubmit.pbs)
		qsub -W depend=afterok:$pvjob -v rad=$rad $curr/pvtrf.pbs
    echo Changing back to working dir: $curr
	cd $curr	
done
