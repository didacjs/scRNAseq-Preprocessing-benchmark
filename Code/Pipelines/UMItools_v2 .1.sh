#!/bin/bash

echo "Starting $(basename "$0")"

mkdir -p output
mkdir -p data
mkdir -p mapping
# Copiar fastqs o descarregarlos
data=$"/home/student/Work/Datasets/5k_mouse_brain_CNIK_3pv3_fastqs"
index=$"/home/student/Work/Refs/STARindex_mmGRCm39/"
t2g=$"/home/student/Work/Refs/t2g/t2g.tsv"
gtf=$"/home/student/Work/Refs/anotacions/Mus_musculus.GRCm39.107.gtf.gz"

# Juntar lanes
setname=$(ls $data/ | sed -E  -n "/[_L00]{4}[0-9]+.*/s/[_L00]{4}[0-9]+.*//p" | head -1)
echo "Setname is $setname" 

# Juntar lanes
numfiles=$(ls $data/ | wc -l)
if (($numfiles == 4))
then 
	echo "No merging needed"
	R1=$(ls $data | grep "R1")
	R1=$(echo $data/$R1)
	R2=$(ls $data | grep "R2")
	R2=$(echo $data/$R2)
else
	if test -f data/READY_"$setname"_R1.fastq.gz
	then
		echo "Lanes already merged"
	else
		echo "Merging lanes"
		mkdir -p data
		cat $data/"$setname"*_R1_*.fastq.gz > data/READY_"$setname"_R1.fastq.gz
		R1="data/READY_"$setname"_R1.fastq.gz"
		cat $data/"$setname"*_R2_*.fastq.gz > data/READY_"$setname"_R2.fastq.gz
		R2=data/READY_"$setname"_R2.fastq.gz
		echo "Lanes merged" 
	fi
fi

# Identificar els BC correctes (whitelist)
if test -f output/whitelist.txt
then
	echo "Whitelist is present"
else
	echo "Running umi_tools whitelist.."
	umi_tools whitelist --stdin $R1 \
        	            --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
        	            --set-cell-number=5000 \
        	            --log2stderr > output/whitelist.txt;
	whitelist=output/whitelist.txt
fi

# Step 3: Extract barcdoes and UMIs and add to read names
if test -f output/EXTRACTED_"$setname"_R1.fastq.gz
then
	echo "Barcodes already extracted"
else
	echo "Running umi_tools extract.."
	umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
		          --stdin $R1 \
		          --stdout data/EXTRACTED_"$setname"_R1.fastq.gz \
		          --read2-in $R2 \
		          --read2-out data/EXTRACTED_"$setname"_R2.fastq.gz \
		          --whitelist output/whitelist.txt;
	R1=data/EXTRACTED_"$setname"_R1.fastq.gz
	R2=data/EXTRACTED_"$setname"_R2.fastq.gz
	
fi

# Step 4: Mapping
if test -d _STARtmp
then
	echo "Mapping had already started but did not finish."
else
	echo "Running STAR mapping.."
	STAR --runThreadN 6 \
	     --genomeDir $index \
	     --readFilesIn $R2 $R1 \
	     --readFilesCommand zcat \
	     --outFilterMultimapNmax 1 \
	     --outSAMtype BAM SortedByCoordinate \
	     --outFileNamePrefix mapping/;
fi

# Step 5: Assign reads to genes
featureCounts -p \
	--countReadPairs \
	-a $gtf \
	-o mapping/gene_assigned.bam \
	-R BAM \
	-T 12 \
	mapping/Aligned.sortedByCoord.out.bam;
	

samtools sort mapping/Aligned.sortedByCoord.out.bam.featureCounts.bam -o mapping/assigned_sorted.bam
samtools index mapping/assigned_sorted.bam 

umi_tools count \
	--per-gene \
	--gene-tag=XT \
	--assigned-status-tag=XS \
	--per-cell \
	-I mapping/assigned_sorted.bam \
	-S output/counts.tsv.gz

echo "End of script."


