#!/bin/bash

cat trackerlist.txt | while read line
do 
	echo ${line#*Data};
	if test -f ~/hpctmp/Data${line#*Data}/VirtualMazeBatchLog.txt
	then
		cat ~/hpctmp/Data${line#*Data}/VirtualMazeBatchLog.txt | tail -n 2; 
		echo;
	else
		echo Log not found;
		echo;
	fi 
done
