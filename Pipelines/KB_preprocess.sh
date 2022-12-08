#!/bin/bash

echo "Starting $(basename "$0")"

source /home/student/Work/Code/scRNAseq-Preprocessing-benchmark/Pipelines/module.sh

mkdir output

# Set dataset and index path
data=$"/home/student/Work/Datasets/Mouse/SRR6835871/10X_P8_15/10X_P8_15_MissingLibrary_1_H3FYJDMXX"
t2g=$"/home/student/Work/Refs/Mouse/T2G/t2g.tsv"
index=$"/home/student/Work/Refs/Mouse/index_KB/index_KB"
tech="10xv2" # Which technolgy was used to produce the dataset

# Find out the name of the dataset
setname=$(ls $data/ | sed -E  -n "/[_L00]{4}[0-9]+.*/s/[_L00]{4}[0-9]+.*//p" | head -1)
echo "Setname is $setname" 

# Merge lanes if needed 
merge_files $data

# Run kb count
echo "Starting count.."
/usr/bin/time -v -o time.txt kb count -i $index -g $t2g -x $tech -o output/ $R1 $R2

# Exit code
exit_code_function "kb count"