#!/bin/bash

txtname='checkmissing_hpc.txt'
cd ~/hpctmp/Data/picasso-misc

for datedir in 2018*; do
    for sessdir in "$datedir"/session0*; do
        if ! [ -f "$sessdir"/vmpv.mat ]; then
            echo "$sessdir"/vmpv.mat >> ~/$txtname
        fi
        if ! [ -f "$sessdir"/binData.hdf ]; then
            echo "$sessdir"/binData.hdf >> ~/$txtname
        fi
        for arraydir in "$sessdir"/array*; do
            for chdir in "$arraydir"/channel*; do
                if [ -d "$chdir"/cell01 ]; then
                    for celldir in "$chdir"/cell*; do
                        if ! [ -f "$celldir"/vmpc.mat ]; then
                            echo "$celldir"/vmpc.mat >> ~/$txtname
                        fi
                        if ! [ -f "$celldir"/vmsv.mat ]; then
                            echo "$celldir"/vmsv.mat >> ~/$txtname
                        fi
                    done
                fi
            done
        done
    done
done

echo "" >> ~/$txtname
