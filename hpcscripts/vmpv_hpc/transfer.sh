#!/bin/bash

curr=$(pwd)
rad=$1
cat sessions.txt | while read line
do
	temp=${line##*picasso-misc}
	echo Transferring files from: ~/hpctmp/Data/picasso-misc$temp
	cd ~/hpctmp/Data/picasso-misc$temp
	scp vmpv.mat hippocampus@cortex.nus.edu.sg:/picasso-misc-folder-link$temp/${rad}vmpv.mat 
done
cd $curr

