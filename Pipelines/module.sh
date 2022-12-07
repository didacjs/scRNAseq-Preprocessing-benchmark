#!/bin/bash

# Takes the name of the function it's evaluating 
exit_code_function () {
	exit_code=$?
	if (($exit_code != 0))
	then
    	echo -e "\nERROR in: $1"
    	exit $exit_code  
	fi
}


# Takes the dataset's directory 
merge_files (){
	if [ ! -d $1 ]
	then
		echo "Directory $1 not found"
		exit
	fi 
	numfiles=$(ls $1/ | wc -l)
	if (($numfiles <= 4))
	then 
		echo "No merging needed"
		R1=$(ls $1 | grep "R1")
		R1=$(echo $1/$R1)
		R2=$(ls $1 | grep "R2")
		R2=$(echo $1/$R2)
		if [ ! -f $R1 ]
		then
			echo "R1 Should be at $R1"
			echo "File $R1 not found"
			exit
		fi
		if [ ! -f $R2 ]
		then
			echo "R2 Should be at $R2"
			echo "File $R2 not found"
			exit
		fi 
	else
		if test -f data/READY_"$setname"_R1.fastq.gz
		then
			echo "Lanes already merged"
		else
			echo "Merging lanes"
			mkdir -p data
			cat $1/"$setname"*_R1_*.fastq.gz > data/READY_"$setname"_R1.fastq.gz
			exit_code_function "cat R1"
			R1=data/READY_"$setname"_R1.fastq.gz
			cat $1/"$setname"*_R2_*.fastq.gz > data/READY_"$setname"_R2.fastq.gz
			exit_code_function "cat R1"
			R2=data/READY_"$setname"_R2.fastq.gz
			echo "Lanes merged" 
		fi
	fi
	
	}
