#!/bin/bash

echo "Starting $(basename "$0")"

source /home/student/Work/Code/scRNAseq-Preprocessing-benchmark/Pipelines/module.sh

mkdir -p output
mkdir -p data
mkdir -p mapping

# Copiar fastqs o descarregarlos
data=$"/home/student/Work/Datasets/Mouse/SRR6835871/10X_P8_15/10X_P8_15_MissingLibrary_1_H3FYJDMXX"
t2g=$"/home/student/Work/Refs/Mouse/T2G/t2g.tsv"
index=$"/home/student/Work/Refs/Mouse/STARindex_mmGRCm39"
gtf=$"/home/student/Work/Refs/Mouse/anotacions/Mus_musculus.GRCm39.107.gtf"

# Get sample name
setname=$(ls $data/ | sed -E  -n "/[_L00]{4}[0-9]+.*/s/[_L00]{4}[0-9]+.*//p" | head -1)
echo "Setname is $setname" 

# Merge lanes if needed 
merge_files $data

# Identify correct barcodes (whitelist)
if test -f output/whitelist.txt
then
	echo "Whitelist is present"
else
	echo "Running umi_tools whitelist.."
	/usr/bin/time -v -o time.txt umi_tools whitelist --stdin $R1 \
							--bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
							--set-cell-number=5000 \
							--log2stderr > output/whitelist.txt;
	whitelist=output/whitelist.txt
fi

exit_code_function "umi_tools whitelist"

# Step 3: Extract barcdoes and UMIs and add to read names
if test -f data/EXTRACTED_"$setname"_R1.fastq.gz
then
	echo "Barcodes already extracted"
else
	echo "Running umi_tools extract.."
	/usr/bin/time -v -a -o  time.txt umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
					--stdin $R1 \
					--stdout data/EXTRACTED_"$setname"_R1.fastq.gz \
					--read2-in $R2 \
					--read2-out data/EXTRACTED_"$setname"_R2.fastq.gz \
					--whitelist output/whitelist.txt \
					--limitBAMsortRAM 11598702704;
	R1=data/EXTRACTED_"$setname"_R1.fastq.gz
	R2=data/EXTRACTED_"$setname"_R2.fastq.gz
	
fi

exit_code_function "umi_tools extract"

# Step 4: Mapping
if test -f mapping/Aligned.out.bam
then
	if test -d mapping/_STARtmp
	then
		echo "Mapping had already started but did not finish."
		exit
	else
		echo "Mapping already performed."
	fi
else
	echo "Running STAR mapping.."
	/usr/bin/time -v -a -o time.txt STAR --runThreadN 6 \
			--genomeDir $index \
			--readFilesIn $R2 $R1 \
			--readFilesCommand zcat \
			--outFilterMultimapNmax 1 \
			--outSAMtype BAM Unsorted \
			--outFileNamePrefix mapping/;
fi

exit_code_function "STAR"

echo "Sorting BAM file"

if test -f mapping/Aligned.Sorted.out.bam
then
	echo "BAM file already sorted"
else
	/usr/bin/time -v -a -o time.txt samtools sort mapping/Aligned.out.bam -o mapping/Aligned.Sorted.out.bam
	exit_code_function "samtools sort"

	/usr/bin/time -v -a -o time.txt samtools index mapping/Aligned.Sorted.out.bam 
	exit_code_function "samtools index"
fi

# Step 5: Assign reads to genes

echo "Assigning reads to genes"
if test -f mapping/Aligned.Sorted.out.bam.featureCounts.bam
then
	echo "Reads already assigned"
else
	/usr/bin/time -v -a -o time.txt featureCounts -p \
			--countReadPairs \
			-R BAM \
			-T 12 \
			-a $gtf \
			-o mapping/gene_assigned.bam \
			mapping/Aligned.Sorted.out.bam;
	exit_code_function "featureCounts"
fi

echo "Sorting assigned BAM file"

if test -f mapping/Assigned.Sorted.out.bam
then
	echo "BAM file already sorted"
else
	/usr/bin/time -v -a -o time.txt samtools sort mapping/Aligned.Sorted.out.bam.featureCounts.bam -o mapping/Assigned.Sorted.out.bam
	exit_code_function "samtools sort"

	/usr/bin/time -v -a -o time.txt samtools index mapping/Assigned.Sorted.out.bam 
	exit_code_function "samtools index"
fi
echo "Running umi_tools count"

/usr/bin/time -v -a -o time.txt umi_tools count \
	--per-gene \
	--gene-tag=XT \
	--assigned-status-tag=XS \
	--per-cell \
	-I mapping/Assigned.Sorted.out.bam \
	-S output/counts.tsv.gz
exit_code_function "umi_tools count"

echo "End of script."



