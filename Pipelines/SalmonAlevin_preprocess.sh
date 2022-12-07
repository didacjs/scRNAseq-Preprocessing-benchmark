#!/bin/bash

echo "Starting $(basename "$0")"

mkdir -p output

# Exit code function
exit_code_function () {
	exit_code=$?
	if (($exit_code != 0))
	then
    	echo -e "\nERROR in: $1"
    	exit $exit_code  
	fi
}

# Copiar fastqs o descarregarlos
data=$"/home/student/Work/Datasets/5k_mouse_brain_CNIK_3pv3_fastqs"
index=$"/home/student/Work/Refs/SALMON/salmon_index"
t2g=$"/home/student/Work/Refs/t2g/t2g.tsv"

# Juntar lanes
setname=$(ls $data/ | sed -E  -n "/[_L00]{4}[0-9]+.*/s/[_L00]{4}[0-9]+.*//p" | head -1)
echo "Setname is $setname" 

# Juntar lanes
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

# Runing salmon alevin

salmon alevin -l ISF \
	-1 $R1 \
	-2 $R2 \
	--chromiumV3 \
	-i $index \
	-p 4 \
	-o output \
	--tgMap $t2g

exit_code_function "salmon alevin"
