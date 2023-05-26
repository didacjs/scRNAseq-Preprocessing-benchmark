#!/bin/bash

echo "Starting $(basename "$0")"

# Script must be executed in the same directory as module.sh
# Contains various functions that are repeated a few times in all pipelines
source "$(dirname -- "${BASH_SOURCE[0]}")/module.sh"

# Usage message
usage="Missing option or argument. Script usage:\n-d: Directory with the FASTQ files.\n-o: Folder to output results.\n-t: Path to the transcript to gene file\n-i: Path to the STAR index directory\n-g: Path to the .gtf anotation file"

# Set up the options and their parameters
while getopts ':d:o:i:t:g:' OPTION; do
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
    
    t)
        t2g=$OPTARG
        if [ ! -f "$t2g" ]
        then 
            echo "No file at provided transcript to gene path"
            echo -e $usage
            exit 2
        fi
        echo "Using t2g from $t2g"
        ;;
	g)
        gtf=$OPTARG
        if [ ! -f "$gtf" ]
        then 
            echo "No file at provided transcript to gene path"
            echo -e $usage
            exit 2
        fi
        echo "Using anotation file from $gtf"
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
if [ -z $t2g ]; then echo "-t not provided"; echo -e $usage; exit 2; fi
if [ -z $gtf ]; then echo "-g not provided"; echo -e $usage; exit 2; fi

exec &> >(tee $output/logfile.txt)

#gtf=$"/home/student/Work/Refs/Mouse/anotacions/Mus_musculus.GRCm39.107.gtf"

# Get sample name
setname=$(ls $data/ | sed -E  -n "/[_L00]{4}[0-9]+.*/s/[_L00]{4}[0-9]+.*//p" | head -1)
echo "Setname is $setname" 

# Merge lanes if needed 
merge_files $data $output/data

# Identify correct barcodes (whitelist)
if test -f $output/whitelist.txt
then
	echo "Whitelist is present"
else
	echo "Running umi_tools whitelist.."
	/usr/bin/time -v -o $output/time.txt umi_tools whitelist --stdin $R1 \
							--bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
							--set-cell-number=6000 \
							-L $output/extract.log \
							--log2stderr > $output/whitelist.txt;
	whitelist=$output/whitelist.txt
fi

exit_code_function "umi_tools whitelist"

# Step 3: Extract barcdoes and UMIs and add to read names
if test -f $output/data/EXTRACTED_"$setname"_R1.fastq.gz
then
	echo "Barcodes already extracted"
	R1=$output/data/EXTRACTED_"$setname"_R1.fastq.gz
	R2=$output/data/EXTRACTED_"$setname"_R2.fastq.gz
else
	echo "R1 is $R1"
	echo "R2 is $R2"
	echo "Running umi_tools extract.."
	/usr/bin/time -v -a -o  $output/time.txt umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
					--stdin $R1 \
					--stdout $output/data/EXTRACTED_"$setname"_R1.fastq.gz \
					--read2-in $R2 \
					--read2-out=$output/data/EXTRACTED_"$setname"_R2.fastq.gz \
					--whitelist=$output/whitelist.txt
	R1=$output/data/EXTRACTED_"$setname"_R1.fastq.gz
	R2=$output/data/EXTRACTED_"$setname"_R2.fastq.gz
	
fi

exit_code_function "umi_tools extract"


echo "R1 is $R1"
echo "R2 is $R2"
# Step 4: Mapping
if test -f $output/mapping/Aligned.out.bam
then
	if test -d $output/mapping/_STARtmp
	then
		echo "Mapping had already started but did not finish."
		exit
	else
		echo "Mapping already performed."
	fi
else
	echo "Running STAR mapping.."
	/usr/bin/time -v -a -o $output/time.txt STAR --runThreadN 6 \
			--genomeDir $index \
			--readFilesIn $R2 \
			--readFilesCommand zcat \
			--outFilterMultimapNmax 1 \
			--outSAMtype BAM Unsorted \
			--outFileNamePrefix $output/mapping/;
fi

exit_code_function "STAR"

# Step 5: Assign reads to genes

echo "Assigning reads to genes"
if test -f $output/mapping/Aligned.out.bam.featureCounts.bam
then
	echo "Reads already assigned"
else
	/usr/bin/time -v -a -o $output/time.txt featureCounts \
			-a $gtf \
			-o $output/mapping/gene_assigned.bam -R BAM $output/mapping/Aligned.out.bam \
			-T 12;
	exit_code_function "featureCounts"
fi

echo "Sorting assigned BAM file"

if test -f $output/mapping/Assigned.Sorted.out.bam
then
	echo "BAM file already sorted"
else
	/usr/bin/time -v -a -o $output/time.txt samtools sort $output/mapping/Aligned.out.bam.featureCounts.bam -o $output/mapping/Assigned.Sorted.out.bam
	exit_code_function "samtools sort"

	/usr/bin/time -v -a -o $output/time.txt samtools index $output/mapping/Assigned.Sorted.out.bam 
	exit_code_function "samtools index"
fi
echo "Running umi_tools count"

/usr/bin/time -v -a -o $output/time.txt umi_tools count \
	--per-gene \
	--gene-tag=XT \
	--assigned-status-tag=XS \
	--per-cell \
	-I $output/mapping/Assigned.Sorted.out.bam \
	-S $output/counts.tsv.gz
exit_code_function "umi_tools count"

echo "End of script."



