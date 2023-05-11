#!/bin/bash

cat ../sessions.txt | while read line
do
	temp=${line##*Data}
	echo Reading from sessions.txt: $line
	dirpath=~/hpctmp/Data$temp
	mkdir -p $dirpath
	echo hpc session directory: $dirpath
	ssh -p 8398 hippocampus@cortex.nus.edu.sg /volume1/Hippocampus/Data/picasso-misc/find_channels.sh $line 
	scp -P 8398 hippocampus@cortex.nus.edu.sg:$line/channels.txt $dirpath
	curr=$(pwd)
	cd $dirpath
		trfjob=$(qsub $curr/mstrfin.pbs)
		msjob=$(qsub -W depend=afterok:$trfjob $curr/mssubmit.pbs)
		qsub -W depend=afterok:$msjob $curr/mstrfout.pbs
	echo Changing back to working dir: $curr
	cd $curr
	ssh -p 8398 hippocampus@cortex.nus.edu.sg rm $line/channels.txt
done
