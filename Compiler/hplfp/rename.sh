#!/bin/sh

#Run line below to rename for all channels in a directory
#cwd=`pwd`; for i in `find . -name "channel*"`; do echo $i; cd $i; $GITHUB_MATLAB/Hippocampus/Compiler/hplfp/rename.sh; cd $cwd; done

while read line; do
	lineArr=($line)
       remoteSize=${lineArr[1]}

      	tmp=${lineArr[0]}
	name=${tmp::${#tmp}-1}
       	fileSize=$(stat -f%z $name)
       	if [ “$fileSize” == “$remoteSize” ] && [ -f $name ]; then 
           deleteName=${tmp::${#tmp}-6}
           deleteName+=.mat
	   rm $deleteName
           mv $name $deleteName
       	fi
done < size.txt

rm size.txt

