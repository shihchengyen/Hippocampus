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
	for file in rplparallel.mat unityfile.mat eyelink.mat umaze.mat
	do
		if [ ! -f $dirpath/$file ]; then
    			scp hippocampus@cortex.nus.edu.sg:$line/$file $dirpath
		fi
	done
	curr=$(pwd)
	cd $dirpath
	echo $dirpath > batch.txt
		rcjob=$(qsub -v rad=$rad $curr/rcsubmit.pbs)
		1pvjob=$(qsub -W depend=afterok:$rcjob -v file_name=1binData.csv $curr/pvsubmit.pbs)
		mvpvjob=$(qsub -W depend=afterok$rcjob -v file_name=mbinData.csv $curr/pvsubmit.pbs)
		qsub -W depend=afterok:$1pvjob -v rad=$rad $curr/rcpvtrf.pbs
		qsub -W depend=afterok:$mpvjob -v rad=$rad $curr/rcpvtrf.pbs
	echo Changing back to working dir: $curr
	cd $curr	
done

