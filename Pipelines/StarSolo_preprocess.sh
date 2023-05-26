#!/bin/bash

echo "Starting $(basename "$0")"

# Script must be executed in the same directory as module.sh
# Contains various functions that are repeated a few times in all pipelines
source "$(dirname -- "${BASH_SOURCE[0]}")/module.sh"

# Usage message
usage="Missing option or argument. Script usage:\n-d: Directory with the FASTQ files.\n-o: Folder to output results.\n-i: Path to the STAR index directory\n-w: Path to the relevant barcode whitelist"

# Set up the options and their parameters
while getopts ':d:o:i:w:' OPTION; do
  case "$OPTION" in
    d)
        data=$OPTARG
        if [ ! -d "$data" ] # Must be provided
        then 
            echo "No data directory at path provided" 
            exit 2
        fi
        echo "Reading FASTQ's from $data"
        ;;
    o)
        output="$OPTARG"
        if [ ! -d "$output" ]
        then
            echo "Creating $output Directory"
            mkdir $output
        fi
        echo "Outputting to $output"
        ;;
    i)
        index=$OPTARG
        if [ ! -d "$index" ]
        then 
            echo "No directory at provided index path"
            echo -e $usage
            exit 2
        fi
        echo "Using index from $index"
        ;;
    
    w)
        whitelist=$OPTARG
        if [ ! -f "$whitelist" ]
        then 
            echo "No file at provided index path"
            echo -e $usage
            exit 2
        fi
        echo "Using index from $index"
        ;;
    ?)
      echo -e $usage
      exit 2
      ;;
  esac
done

# All options are mandatory and the path indicated must exist. 
# The exceptions are -o, in which cas we will create the directory and -T which has a default.
if [ -z $data ]; then echo "-d not provided"; echo -e $usage; exit 2; fi
if [ -z $index ]; then echo "-i not provided"; echo -e $usage; exit 2; fi
if [ -z $output ]; then echo "-o not provided"; echo -e $usage; exit 2; fi
if [ -z $whitelist ]; then echo "-w not provided"; echo -e $usage; exit 2; fi

exec &> >(tee $output/logfile.txt)

# Indicar ruta del dataset, el transcript to gene i l'index kallisto de l'especie
# data=$"/home/student/Work/Datasets/Mouse/5k_mouse_brain_CNIK_3pv3_fastqs/"
# whitelist=$"/home/student/Work/Refs/Mouse/Whitelists/737K-august-2016.txt"
# index=$"/home/student/Work/Refs/Mouse/STARindex_mmGRCm39"

# Nom del dataset
setname=$(ls $data/ | sed -E  -n "/[_L00]{4}[0-9]+.*/s/[_L00]{4}[0-9]+.*//p" | head -1)
echo "Setname is $setname" 

# Merge lanes if needed 
merge_files $data $output/data

# Unzip
echo "Gunziping $R1"
gunzip $R1
echo "Gunziping $R2"
gunzip $R2
R1=$output/data/READY_"$setname"_R1.fastq
R2=$output/data/READY_"$setname"_R2.fastq

exit_code_function "gunzip"

# Run STAR with --soloType set to CB_UMI_Simple
/usr/bin/time -v -a -o $output/time.txt STAR --genomeDir $index \
     --readFilesIn $R2 $R1 \
     --soloType CB_UMI_Simple \
     --soloCBwhitelist  $whitelist \
     --outFileNamePrefix $output/ \
     --runThreadN 6

# Exit code
exit_code_function "STAR"