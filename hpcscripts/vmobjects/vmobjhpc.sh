#!/bin/bash

cat ../cells.txt | while read line
do
	celldir=${line##*Data}
	echo Reading from cells.txt: $line
	sessdir=~/hpctmp/Data${celldir%%/array*}
	mkdir -p $sessdir
	echo hpc session directory: $sessdir
	for file in 1vmpv.mat rplparallel.mat unityfile.mat umaze.mat
	do
        if [ ! -f $sessdir/$file ]; then
            scp hippocampus@cortex.nus.edu.sg:${line%%/array*}/$file $sessdir
        fi
	done
	dirpath=~/hpctmp/Data$celldir
	mkdir -p $dirpath
	echo hpc cell directory: $dirpath
    if [ ! -f $dirpath/spiketrain.mat ]; then
        scp hippocampus@cortex.nus.edu.sg:$line/spiketrain.mat $dirpath
    fi
	for file in rplparallel.mat unityfile.mat umaze.mat
	do
		if [ ! -f $dirpath/$file ]; then
			ln -s $sessdir/$file $dirpath/$file
		fi
	done
	curr=$(pwd)
    cd $dirpath
        pcjob=$(qsub $curr/vmpsubmit.pbs)
        svjob=$(qsub $curr/vmssubmit.pbs)
        hdjob=$(qsub $curr/vmhsubmit.pbs)
        qsub -W depend=afterok:$pcjob:$svjob:$hdjob -v dir=$dirpath $curr/vmobjtrf.pbs
    echo Changing back to working dir: $curr
	cd $curr
    echo $dirpath > vmobjbatch.txt
done

