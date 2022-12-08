#!/bin/bash

echo "Starting $(basename "$0")"

source module.sh

mkdir -p output
mkdir -p data
mkdir -p mapping

# Copiar fastqs o descarregarlos
data=$"/home/student/Work/Datasets/5k_mouse_brain_CNIK_3pv3_fastqs"
index=$"/home/student/Work/Refs/STARindex_mmGRCm39/"
t2g=$"/home/student/Work/Refs/t2g/t2g.tsv"
gtf=$"/home/student/Work/Refs/anotacions/Mus_musculus.GRCm39.107.gtf.gz"

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
if test -f output/EXTRACTED_"$setname"_R1.fastq.gz
then
	echo "Barcodes already extracted"
else
	echo "Running umi_tools extract.."
	/usr/bin/time -v -o -a time.txt umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
					--stdin $R1 \
					--stdout data/EXTRACTED_"$setname"_R1.fastq.gz \
					--read2-in $R2 \
					--read2-out data/EXTRACTED_"$setname"_R2.fastq.gz \
					--whitelist output/whitelist.txt;
	R1=data/EXTRACTED_"$setname"_R1.fastq.gz
	R2=data/EXTRACTED_"$setname"_R2.fastq.gz
	
fi

exit_code_function "umi_tools extract"

# Step 4: Mapping
if test -d _STARtmp
then
	echo "Mapping had already started but did not finish."
else
	echo "Running STAR mapping.."
	/usr/bin/time -v -o -a time.txt STAR --runThreadN 6 \
			--genomeDir $index \
			--readFilesIn $R2 $R1 \
			--readFilesCommand zcat \
			--outFilterMultimapNmax 1 \
			--outSAMtype BAM SortedByCoordinate \
			--outFileNamePrefix mapping/;
fi

exit_code_function "STAR"

# Step 5: Assign reads to genes
/usr/bin/time -v -o -a time.txt featureCounts -p \
		--countReadPairs \
		-a $gtf \
		-o mapping/gene_assigned.bam \
		-R BAM \
		-T 12 \
		mapping/Aligned.sortedByCoord.out.bam;
exit_code_function "featureCounts"

/usr/bin/time -v -o -a time.txt samtools sort mapping/Aligned.sortedByCoord.out.bam.featureCounts.bam -o mapping/assigned_sorted.bam
exit_code_function "samtools sort"

/usr/bin/time -v -o -a time.txt samtools index mapping/assigned_sorted.bam 
exit_code_function "samtools index"

/usr/bin/time -v -o -a time.txt umi_tools count \
	--per-gene \
	--gene-tag=XT \
	--assigned-status-tag=XS \
	--per-cell \
	-I mapping/assigned_sorted.bam \
	-S output/counts.tsv.gz
exit_code_function "umi_tools count"

echo "End of script."


