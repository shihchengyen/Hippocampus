#!/bin/bash

rad=$1
rad_str=$(echo "$rad" | awk '{gsub(/\./,"-",$1); print $1}')
prefix="r${rad_str}"
vmpv_name="${prefix}vmpv.mat"

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
    mvms=$(qsub -v prefix="$prefix" $curr/vmssubmit.pbs)
    singlevms=$(qsub -v prefix="" $curr/vmssubmit.pbs)
    #echo "qsub -v prefix="$prefix" -W dependon=afterok:"$mvms" "$curr/plotvms.pbs""
    mplot=$(qsub -v prefix="$prefix" -W depend=afterok:"$mvms" "$curr/plotvms.pbs")
    singleplot=$(qsub -v prefix="" -W depend=afterok:"$singlevms" "$curr/plotvms.pbs")
    echo Changing back to working dir: $curr
    cd $curr
    echo ~/hpctmp/Data$celldir > vmsvbatch.txt
done
