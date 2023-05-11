#!/bin/bash

cat ../cells.txt | while read line
do
	celldir=${line##*Data}
	echo Reading from cells.txt: $line
    sessdir=~/hpctmp/Data${celldir%%/array*}
    mkdir -p $sessdir
    echo hpc session directory: $sessdir
    if [ ! -f $sessdir/1vmpv.mat ]; then
        scp -P 8398 hippocampus@cortex.nus.edu.sg:${line%%/array*}/1vmpv.mat $sessdir
    fi
	if [ ! -f $sessdir/rplparallel.mat ]; then
        scp -P 8398 hippocampus@cortex.nus.edu.sg:${line%%/array*}/rplparallel.mat $sessdir
    fi
	if [ ! -f $sessdir/unityfile.mat ]; then
        scp -P 8398 hippocampus@cortex.nus.edu.sg:${line%%/array*}/unityfile.mat $sessdir
    fi
	if [ ! -f $sessdir/umaze.mat ]; then
        scp -P 8398 hippocampus@cortex.nus.edu.sg:${line%%/array*}/umaze.mat $sessdir
    fi
	dirpath=~/hpctmp/Data$celldir
    mkdir -p $dirpath
	echo hpc cell directory: $dirpath
    if [ ! -f $dirpath/spiketrain.mat ]; then
        scp -P 8398 hippocampus@cortex.nus.edu.sg:$line/spiketrain.mat $dirpath
    fi
	if [ ! -f $dirpath/rplparallel.mat ]; then
		ln -s $sessdir/rplparallel.mat $dirpath/rplparallel.mat
	fi
	if [ ! -f $dirpath/unityfile.mat ]; then
		ln -s $sessdir/unityfile.mat $dirpath/unityfile.mat
	fi
	if [ ! -f $dirpath/umaze.mat ]; then
		ln -s $sessdir/umaze.mat $dirpath/umaze.mat
	fi
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
