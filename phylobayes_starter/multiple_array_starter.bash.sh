#!/bin/bash
#
# Usage:
# 	sh ./multiple_array_starter.bash.sh DIRECTORY EXTENSION
#
# Start sbatch array script for each subdirectory found in DIRECTORY
#

if [ "$#" -ne 2 ]; then
	echo "Invalid number of arguments."
	echo "Usage:"
	echo "sh ./multiple_array_starter.bash.sh DIRECTORY EXTENSION"
	exit 1
fi

for directory in `ls -d $1`; do
    sbatch start_phylobayes.sbatch_array.sh $directory $2
    sleep 10
done

