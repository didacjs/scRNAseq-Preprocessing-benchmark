#!/bin/bash

echo "Starting $(basename "$0")"

source /home/student/Work/Code/scRNAseq-Preprocessing-benchmark/Pipelines/module.sh

mkdir output

# Copiar fastqs o descarregarlos
data=$"/home/student/Work/Datasets/Mouse/SRR6835871/10X_P8_15/10X_P8_15_MissingLibrary_1_H3FYJDMXX"
index=$"/home/student/Work/Refs/Mouse/indexes/SA/salmon_index"
t2g=$"/home/student/Work/Refs/Mouse/t2g_gencode_mouse32.tsv"
whitelist="/home/student/Work/Code/scRNAseq-Preprocessing-benchmark/Refs/TM_BC.tsv" # Use BC idenwhitelist
# Get sample name
setname=$(ls $data/ | sed -E  -n "/[_L00]{4}[0-9]+.*/s/[_L00]{4}[0-9]+.*//p" | head -1)
echo "Setname is $setname" 

# Merge lanes if needed 
merge_files $data

# Runing salmon alevin

/usr/bin/time -v -o time.txt salmon alevin \
		-l ISF \
		-1 $R1 \
		-2 $R2 \
		--chromium \
		-i $index \
		-p 4 \
		-o output \
		--tgMap $t2g \
		#--whitelist $whitelist \ 
exit_code_function "salmon alevin"
