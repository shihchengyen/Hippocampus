#!/bin/bash

txtname='moveerrors22june.txt'
movedir='/Volumes/Hippocampus/Data/picasso-misc'
cd /Volumes/Drive1/hpc_processing/22june_batch/Data/picasso-misc

for datedir in 201806*; do
    for sessdir in "$datedir"/session0*; do
        if [ -f "$sessdir"/vmpv.mat ]; then
            cp "$sessdir"/vmpv.mat "$movedir"/"$sessdir"/1vmpv.mat
        else
            echo "$sessdir"/vmpv.mat >> ~/Documents/$txtname
        fi
        if [ -f "$sessdir"/binData.hdf ]; then
            cp "$sessdir"/binData.hdf "$movedir"/"$sessdir"/1binData.hdf
        else
            echo "$sessdir"/binData.hdf >> ~/Documents/$txtname
        fi
        for arraydir in "$sessdir"/array*; do
            for chdir in "$arraydir"/channel*; do
                if [ -d "$chdir"/cell01 ]; then
                    for celldir in "$chdir"/cell*; do
                        if [ -f "$celldir"/vmpc.mat ]; then
                            mkdir -p "$movedir"/"$celldir"/FiltVel/1px
                            cp "$celldir"/vmpc.mat "$movedir"/"$celldir"/FiltVel/1px/vmpc.mat
                        else
                            echo "$celldir"/vmpc.mat >> ~/Documents/$txtname
                        fi
                        if [ -f "$celldir"/vmsv.mat ]; then
                            mkdir -p "$movedir"/"$celldir"/FiltVel/1px
                            cp "$celldir"/vmsv.mat "$movedir"/"$celldir"/FiltVel/1px/vmsv.mat
                        else
                            echo "$celldir"/vmsv.mat >> ~/Documents/$txtname
                        fi
                    done
                fi
            done
        done
    done
    echo $datedir done
done
