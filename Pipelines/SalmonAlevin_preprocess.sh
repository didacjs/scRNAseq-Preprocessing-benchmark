#!/bin/bash

echo "Starting $(basename "$0")"

source module.sh

mkdir output

# Copiar fastqs o descarregarlos
data=$"/home/student/Work/Datasets/5k_mouse_brain_CNIK_3pv3_fastqs"
index=$"/home/student/Work/Refs/SALMON/salmon_index"
t2g=$"/home/student/Work/Refs/t2g/t2g.tsv"

# Get sample name
setname=$(ls $data/ | sed -E  -n "/[_L00]{4}[0-9]+.*/s/[_L00]{4}[0-9]+.*//p" | head -1)
echo "Setname is $setname" 

# Merge lanes if needed 
merge_files $data

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
