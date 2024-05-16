#!/bin/bash
# Check if the radii.txt file exists
if [ ! -f radii.txt ]; then
  echo "Error: radii.txt file not found!"
  exit 1
fi

# Loop through each line in the rad_values.txt file
while IFS= read -r rad || [ -n "$rad" ]; do
  echo "Invoking rchpc.sh with rad: $rad"
  ./vmshpc.sh "$rad"
done < radii.txt

