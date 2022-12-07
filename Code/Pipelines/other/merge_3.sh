#!/bin/bash

listSamples=()

setname=$(ls data/ | grep -E -o '^[^\w L00]+' | head -1)
echo "Set name is $setname"

# loop to add the unique samples to an array
for i in data/*.fastq.gz;

do
    cat data/"$setname"*_R1_* > data/MERGED"$setname"_R1.fastq.gz
    cat data/"$setname"*_R2_* > data/MERGED"$setname"_R2.fastq.gz
    cat data/"$setname"*_I1_* > data/MERGED"$setname"_I1.fastq.gz
    cat data/"$setname"*_I2_* > data/MERGED"$setname"_I2.fastq.gz
    
done
echo "UUUUH"
