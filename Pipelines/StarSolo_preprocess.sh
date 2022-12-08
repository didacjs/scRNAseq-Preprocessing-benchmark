#!/bin/bash

echo "Starting $(basename "$0")"

source /home/student/Work/Code/scRNAseq-Preprocessing-benchmark/Pipelines/module.sh

mkdir output

# Indicar ruta del dataset, el transcript to gene i l'index kallisto de l'especie
data=$"/home/student/Work/Datasets/Mouse/5k_mouse_brain_CNIK_3pv3_fastqs/"
whitelist=$"/home/student/Work/Refs/Mouse/Whitelists/737K-august-2016.txt"
index=$"/home/student/Work/Refs/Mouse/STARindex_mmGRCm39"

# Nom del dataset
setname=$(ls $data/ | sed -E  -n "/[_L00]{4}[0-9]+.*/s/[_L00]{4}[0-9]+.*//p" | head -1)
echo "Setname is $setname" 

# Merge lanes if needed 
merge_files $data

# Unzip
gunzip $R1
exit_code_function "gunzip"
gunzip $R2
exit_code_function "gunzip"
R1=data/READY_"$setname"_R1.fastq
R2=data/READY_"$setname"_R2.fastq

exit_code_function "gunzip"

# Run STAR with --soloType set to CB_UMI_Simple
STAR --genomeDir $index \
     --readFilesIn $R2 $R1 \
     --soloType CB_UMI_Simple \
     --soloCBwhitelist  $whitelist \
     --soloUMIlen 12 \
     --runThreadN 4

# Exit code
exit_code_function "STAR"