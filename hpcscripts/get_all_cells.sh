#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <directory_name> <ssh_address>"
    exit 1
fi

# Assign arguments to variables
dir_name="$1"
ssh_address="$2"

# SSH and execute the command
ssh "$ssh_address" << EOF | while IFS= read -r line; do
    ls -1 -d /${dir_name}/array*/channel*/cell* -f
EOF
    echo "$line"
done

