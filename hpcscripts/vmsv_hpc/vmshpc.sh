#!/bin/bash

cat cells.txt | while read line
do
	celldir=${line##*Data}
	echo Reading from cells.txt: $line
    mkdir -p ~/hpctmp/Data$celldir
	echo hpc cell directory: ~/hpctmp/Data$celldir
    spiketrainFILE=~/hpctmp/Data"$celldir"/spiketrain.mat
    if [ ! -f $spiketrainFILE ]; then
	echo "Copying spiketrain from hippocampus"
        scp hippocampus@cortex.nus.edu.sg:$line/spiketrain.mat ~/hpctmp/Data$celldir
    fi
    
    sessdir=${celldir%%/array*}
    hippsessdir=${line%%/array*}
    mkdir -p ~/hpctmp/Data$sessdir
    echo hippocampus session directory: $hippsessdir
    echo hpc session directory: ~/hpctmp/Data$sessdir
    vmpvFILE=~/hpctmp/Data"$sessdir"/vmpv.mat
    if [ ! -f $vmpvFILE ]; then
        scp hippocampus@cortex.nus.edu.sg:$hippsessdir/1vmpv.mat ~/hpctmp/Data$sessdir/vmpv.mat
    fi
    
    curr=$(pwd)
    cd ~/hpctmp/Data$celldir
        qsub $curr/vmssubmit.txt
    echo Changing back to working dir: $curr
    cd $curr
    echo ~/hpctmp/Data$celldir > vmsvbatch.txt
done
