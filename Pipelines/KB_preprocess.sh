#!/bin/bash

echo "Starting $(basename "$0")"

mkdir output

# Exit code function
exit_code_function () {
	exit_code=$?
	if (($exit_code != 0))
	then
    	echo -e "\nERROR in: $1"
    	exit $exit_code  
	fi
}

# Set dataset and index path
data=$"/home/student/Work/Datasets/Human/GSM3295024/"
t2g=$"/home/student/Work/Refs/Human/T2G/t2g.tsv"
index=$"/home/student/Work/Refs/Human/index_KB/index_KB"

# Which technolgy was used to produce the dataset
tech="10xv3"

# Find out the name of the dataset
setname=$(ls $data/ | sed -E  -n "/[_L00]{4}[0-9]+.*/s/[_L00]{4}[0-9]+.*//p" | head -1)
echo "Setname is $setname" 

# Merge lanes if needed 
numfiles=$(ls $data/ | wc -l)
if (($numfiles <= 4))
then 
	echo "No merging needed"
	R1=$(ls $data | grep "R1")
	R1=$(echo $data/$R1)
	R2=$(ls $data | grep "R2")
	R2=$(echo $data/$R2)
else
	if test -f data/READY_"$setname"_R1.fastq.gz
	then
		echo "Lanes already merged"
	else
		echo "Merging lanes"
		mkdir -p data
		cat $data/"$setname"*_R1_*.fastq.gz > data/READY_"$setname"_R1.fastq.gz
		R1=data/READY_"$setname"_R1.fastq.gz
		cat $data/"$setname"*_R2_*.fastq.gz > data/READY_"$setname"_R2.fastq.gz
		R2=data/READY_"$setname"_R2.fastq.gz
		echo "Lanes merged" 
	fi
fi



# Run kb count
echo "Starting count.."
kb count -i $index -g $t2g -x $tech -o output/ $R1 $R2

# Exit code
exit_code_function "kb count"