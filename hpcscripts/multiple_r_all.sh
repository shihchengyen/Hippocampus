#!/bin/bash
# Check if the radii.txt file exists
if [ ! -f radii.txt ]; then
  echo "Error: radii.txt file not found!"
  exit 1
fi

# Loop through each line in the rad_values.txt file
while IFS= read -r rad || [ -n "$rad" ]; do
  echo "Invoking working.sh with rad: $rad"
  ./working.sh "$rad"
done < radii.txt

