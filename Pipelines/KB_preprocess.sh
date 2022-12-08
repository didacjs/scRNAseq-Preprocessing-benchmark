#!/bin/bash

echo "Starting $(basename "$0")"

source module.sh

mkdir output

# Set dataset and index path
data=$"/media/student/Storage/Datasets/1k_mouse_kidney_CNIK_3pv3_fastqs"
t2g=$"/home/student/Work/Refs/Mouse/T2G/t2g.tsv"
index=$"/home/student/Work/Refs/Mouse/index_KB/index_KB"

# Which technolgy was used to produce the dataset
tech="10xv3"

# Find out the name of the dataset
setname=$(ls $data/ | sed -E  -n "/[_L00]{4}[0-9]+.*/s/[_L00]{4}[0-9]+.*//p" | head -1)
echo "Setname is $setname" 

# Merge lanes if needed 
merge_files $data

# Run kb count
echo "Starting count.."
kb count -i $index -g $t2g -x $tech -o output/ $R1 $R2

# Exit code
exit_code_function "kb count"