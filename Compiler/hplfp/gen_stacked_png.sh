#!/bin/sh

# Path
cd /volume1/Hippocampus/Data ;

# Find array directories with png files generated after certain date
find . -mindepth 4 -maxdepth 4 -type d -name "array*" | grep -v "analog" | cut -d "/" -f 2- | sort > array.txt ;

# Traverse channel directories and create png.tar which contains data of lfp.png and hp.png
for i in $(cat array.txt) ; # Array directory
do
  cd $i ;
  echo $i ;
  find . -mindepth 2 -maxdepth 2 -type f -name "*.png" | cut -d "/" -f 2- | grep -v "unit" | xargs tar -cvf png.tar ;
  cd $cwd ;
done ;
