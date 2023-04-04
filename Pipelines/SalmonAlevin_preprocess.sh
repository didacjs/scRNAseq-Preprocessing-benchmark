#!/bin/bash

echo "Starting $(basename "$0")"

# Script must be executed in the same directory as module.sh
# Contains various functions that are repeated a few times in all pipelines
source "$(dirname -- "${BASH_SOURCE[0]}")/module.sh"

# Usage message
usage="Missing option or argument. Script usage:\n-d: Directory with the FASTQ files.\n-o: Folder to output results.\n-t: Path to the transcript to gene file, with gencode nomenclature of version matching index\n-i: Path to the Salmon Alevin index directory\n-w: Path to the relevant barcode whitelist"

# Set up the options and their parameters
while getopts ':d:o:t:i:w:' OPTION; do
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
    t)
        t2g=$OPTARG
        if [ ! -f "$t2g" ]
        then 
            echo "No file at provided transcript to gene path"
            echo -e $usage
            exit 2
        fi
        echo "Using transcript to gene file from $t2g"
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
if [ -z $t2g ]; then echo "-t not provided"; echo -e $usage; exit 2; fi
if [ -z $index ]; then echo "-i not provided"; echo -e $usage; exit 2; fi
if [ -z $output ]; then echo "-o not provided"; echo -e $usage; exit 2; fi
if [ -z $whitelist ]; then echo "-w not provided"; echo -e $usage; exit 2; fi

exec &> >(tee $output/logfile.txt)

# Get sample name
setname=$(ls $data/ | sed -E  -n "/[_L00]{4}[0-9]+.*/s/[_L00]{4}[0-9]+.*//p" | head -1)
echo "Setname is $setname" 

# Merge lanes if needed 
merge_files $data $output/data

# Runing salmon alevin

/usr/bin/time -v -o $output/time.txt salmon alevin \
		-l ISF \
		-1 $R1 \
		-2 $R2 \
		--chromium \
		-i $index \
		-p 4 \
		-o $output/output \
		--tgMap $t2g \
		#--whitelist $whitelist \ 
exit_code_function "salmon alevin"
