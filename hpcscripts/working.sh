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
    #echo "qsub -v rad="$rad" "$curr/raycasting/rcsubmit.pbs" -N "${line}_${multi_output_file}_raycast""
    rcjob=$(qsub -v rad="$rad" -N "${rad}_raycast" "$curr/raycasting/rcsubmit.pbs" )
    echo $rcjob
    qsub -W depend=afterok:"$rcjob" -v rad="$rad" "$curr/raycasting/rctrf.pbs"
    echo "curr : $curr"
    binnerjob=$(qsub -W depend=afterok:"$rcjob"  -N "${rad}_bin" -v curr="$curr",single_path="unityfile_eyelink.csv",single_save_path="1binData.csv",multi_path="$multi_output_file",multi_save_path="$multi_save_file"  "$curr/raycasting/binning.pbs")

    echo "Changing back to working dir: $curr"
    cd "$curr"

    echo "Setting up vmpv"
    pvjob1=$(qsub -v file_name=1binData.csv -W depend=afterok:"$binnerjob"  -N "1_pv" "$curr/vmpv_hpc/pvsubmit.pbs")
    qsub -W depend=afterok:"$pvjob1" -v prefix="$prefix" "$curr/vmpv_hpc/pvtrf.pbs"
    echo "Setting vmpv ${multi_output_file}"
    pvjobm=$(qsub -v file_name=$multi_save_file -W depend=afterok:"$binnerjob" -N "${rad}_pv" "$curr/vmpv_hpc/pvsubmit.pbs")
    qsub -W depend=afterok:"$pvjobm" -v prefix="$prefix" "$curr/vmpv_hpc/pvtrf.pbs"
    echo Changing back to working dir: "$curr"
    cd "$curr"


    rad_str=$(echo "$rad" | awk '{gsub(/\./,"-",$1); print $1}')
    prefix="r${rad_str}"
    vmpv_name="${prefix}vmpv.mat"
    
    echo "Doing vmsv"
    $curr/get_all_cells.sh "$line" "hippocampus@cortex.nus.edu.sg" | while IFS= read -r cell; do
        celldir=${cell##*Data}
        mkdir -p ~/hpctmp/Data"$celldir"
        
        spiketrainFILE=~hpctmp/Data"$celldir"/spiketrain.mat
        if [ ! -f $spiketrainFILE ]; then
	    echo "Copying spiketrain from hippocampus"
            scp hippocampus@cortex.nus.edu.sg:$cell/spiketrain.mat ~/hpctmp/Data$celldir
        fi
        cd ~/hpctmp/Data$celldir
        echo "${cell}_${prefix}vms"
        echo "PV job : $pvjobm , $pvjob1" 
        mvms=$(qsub -v prefix="$prefix" -N "${prefix}vms" -W depend=afterok:"$pvjobm" $curr/vmsv_hpc/vmssubmit.pbs) 
        singlevms=$(qsub -v prefix="" -N "1vms" -W depend=afterok:"$pvjob1" $curr/vmsv_hpc/vmssubmit.pbs) 
        qsub -W depend=afterok:"$mvms" -v prefix="$prefix" "$curr/vmsv_hpc/svtrf.pbs"
        qsub -W depend=afterok:"$singlevms" -v prefix="$prefix" "$curr/vmsv_hpc/svtrf.pbs"
	### does not seem to be working for some reason
	### mplot=$(qsub -v prefix="$prefix" -N "${prefix}plot" -W depend=afterok:"$mvms" "$curr/vmsv_hpc/plotvms.pbs")
        ### singleplot=$(qsub -v prefix="" -N "singleplot" -W depend=afterok:"$singlevms" "$curr/vmsv_hpc/plotvms.pbs")
        echo Changing back to working dir: $curr
        cd $curr
     done

done < sessions.txt

