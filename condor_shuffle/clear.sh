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
		rm *.0.*
		rm vmplacecell*
		cd $ori
	fi
done < "$file"
