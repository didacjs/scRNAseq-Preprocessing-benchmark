#!/bin/bash

echo "Starting $(basename "$0")"

source /home/student/Work/Code/scRNAseq-Preprocessing-benchmark/Pipelines/module.sh

mkdir output

# Copiar fastqs o descarregarlos
data=$"/home/student/Work/Datasets/Mouse/SRR6835871/10X_P8_15/10X_P8_15_MissingLibrary_1_H3FYJDMXX"
index=$"/home/student/Work/Refs/Mouse/SALMON/salmon_index"
t2g=$"/home/student/Work/Refs/Mouse/T2G/t2g.tsv"

# Get sample name
setname=$(ls $data/ | sed -E  -n "/[_L00]{4}[0-9]+.*/s/[_L00]{4}[0-9]+.*//p" | head -1)
echo "Setname is $setname" 

# Merge lanes if needed 
merge_files $data

# Runing salmon alevin

/usr/bin/time -v -o time.txt salmon alevin -l ISF \
		-1 $R1 \
		-2 $R2 \
		--chromium \
		-i $index \
		-p 4 \
		-o output \
		--tgMap $t2g \
		--expectCells 12000 # Cells found by tabula muris preprocess

exit_code_function "salmon alevin"
