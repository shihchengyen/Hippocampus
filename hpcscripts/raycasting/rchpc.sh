#!/bin/bash
#set -x

rad=$1
rad_str=$(echo "$rad" | awk '{gsub(/\./,"-",$1); print $1}')
multi_output_file="unityfile_eyelink_r${rad_str}.csv"
multi_save_file="r${rad_str}binData.csv"
while IFS= read -r line || [ -n "$line" ]; do
    temp="${line##*Data}/session01"
    echo "Reading from sessions.txt: $line"
    line+="/session01"
    dirpath=~/hpctmp/Data$temp
    mkdir -p "$dirpath"
    echo "---"
    echo "hpc session directory: $dirpath"
    for file in unityfile.mat eyelink.mat; do
        if [ ! -f "$dirpath/$file" ]; then
            scp "hippocampus@cortex.nus.edu.sg:$line/$file" "$dirpath"
        fi
    done

    curr=$(pwd)
    cd "$dirpath"
    echo "$dirpath" > batch.txt
    rcjob=$(qsub -v rad="$rad" "$curr/rcsubmit.pbs")
    qsub -W depend=afterok:"$rcjob" -v rad="$rad" "$curr/rctrf.pbs"
    echo "curr : $curr"
    binnerjob=$(qsub -W depend=afterok:"$rcjob" -v curr="$curr",single_path="unityfile_eyelink.csv",single_save_path="binData.csv",multi_path="$multi_output_file",multi_save_path="$multi_save_file"  "$curr/binning.pbs")

    echo "Changing back to working dir: $curr"
    cd "$curr"
done < batch.txt

