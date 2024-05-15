#!/bin/bash

rad=$1
rad_str=$(echo "$rad" | awk '{gsub(/\./,"-",$1); print $1}')
bin_file_name="r${rad_str}binData.csv"

cat ../sessions.txt | while read line
do
        temp="${line##*Data}/session01"
        echo Reading from sessions.txt: $line
        line+="/session01"
        dirpath=~/hpctmp/Data$temp
        mkdir -p "$dirpath"
        echo hpc session directory: "$dirpath"
        for file in mbinData.csv 1binData.csv rplparallel.mat umaze.mat
        do
                if [ ! -f "$dirpath/$file" ]; then
                        scp hippocampus@cortex.nus.edu.sg:"$line/$file" "$dirpath"
                fi
        done
        curr=$(pwd)
        cd "$dirpath"
        echo "$dirpath" > batch.txt
        pvjob1=$(qsub -v file_name=binData.csv "$curr/pvsubmit.pbs")
        qsub -W depend=afterok:"$pvjob1" -v rad="$rad" "$curr/pvtrf.pbs"
	echo "Getting ${bin_file_name}"
        pvjobm=$(qsub -v file_name=$bin_file_name "$curr/pvsubmit.pbs")
        qsub -W depend=afterok:"$pvjobm" -v rad="$rad" "$curr/pvtrf.pbs"
        echo Changing back to working dir: "$curr"
        cd "$curr"
done

