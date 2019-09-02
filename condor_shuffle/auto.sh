#!/bin/sh

ori=$(pwd)

file="cell_list.txt"
while read line
do
	size=${#line}
	if [ $size -gt 10 ]
	then
		echo $line
		cd $line
		condor_submit /Volumes/User/ngkianwei/Desktop/condor_shuffle/shuffle_submit.txt	
		cd $ori
	fi
done < "$file"
