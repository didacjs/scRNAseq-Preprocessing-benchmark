#!/bin/bash

echo "Starting $(basename "$0")"

# Script must be executed in the same directory as module.sh
# Contains various functions that are repeated a few times in all pipelines
source "$(dirname -- "${BASH_SOURCE[0]}")/module.sh"

# Usage message
usage="Missing option or argument. Script usage:\n-d: Directory with the FASTQ files.\n-o: Folder to output results.\n-t: Path to the transcript to gene file\n-i: Path to the Kallisto Bustools index file\n[-T: scRNA-Seq technology, a name from kb --list. Default is 10xv2]"

# Set up the options and their parameters
while getopts ':d:o:t:i:T:' OPTION; do
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
        if [ ! -f "$index" ]
        then 
            echo "No file at provided index path"
            echo -e $usage
            exit 2
        fi
        echo "Using index from $index"
        ;;
    
    T)
        tech=$OPTARG
        list=$(kb --list | awk '{print $1}' | tail -18)
        match=$(echo $list | grep -ow $tech)
        if [ ! $tech = $match ]
        then 
            echo "Tech provided not supported"
            echo -e $usage
            exit 2
        fi
        echo "Using kallisto bustools for $tech"
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
if [ -z $tech ]; then echo "Using Illumina 10xv2 parameters"; tech=10xv2; fi

exec &> >(tee $output/logfile.txt)

# Find out the name of the dataset
setname=$(ls $data/ | sed -E  -n "/[_L00]{4}[0-9]+.*/s/[_L00]{4}[0-9]+.*//p" | head -1)
echo "Setname is $setname" 

# Merge lanes if needed, save merged FASTQs at $output/data
merge_files $data $output/data

echo $R1
echo $R2
# Run kb count
echo "Starting count.."
/usr/bin/time -v -o $output/time.txt kb count -i $index -g $t2g -x $tech -o $output/Counts $R1 $R2 --cellranger

# Exit code
exit_code_function "kb count"