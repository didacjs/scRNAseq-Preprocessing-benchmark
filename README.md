# Benchmark and downstream comprative of four scRNA-Seq preprocessing pipelines (W.I.P.)

## General description

Single cell transcriptomic techniques (scRNA-Seq) examine the gene expression level in individual cells, simultaneously measuring concentration of messenger RNA transcripts (mRNA) belonging to hundreds to thousands of genes. The analysis of transcriptomics at this  resolution makes possible the decryption of heterogeneous cell populations, the study of cellular life cycles and the modeling of transcriptional dynamics.

The raw output of these techniques are sequencing files (FASTQ) that contain millions of sequences (reads) of transcripts. Each read contains a nucleotide sequence encoding a barcode and a Unique Molecular Identifier (UMI). The barcode associates the read with a cell and the UMI with an RNA molecule. These files cannot be analyzed directly. As a first step, it is necessary to preprocess the raw data to produce a counts matrix, which contains the number of  transcripts for each gene in each cell. With the counts matrix, we can then study the transcriptomic profile of our cell sample. 

The preprocessing step can be computationally intensive depending on the size of the sample, and so a variety of specialized tools to carry and optimize the task have been developed. Despite that, there is still no procedure considered as a standard. Consequently, it is interesting to compare the different available tools to find the one that best suits our needs, computational resources and specifications of the experiment.

I have designed four scRNA-seq preprocessing pipelines in bash, each starting from one different tool, STAR/UMI-Tools, STARsolo, Kallisto|Bustools and Salmon Alevin. These pipelines take FASTQ files as input and generate a counts matrix, along with a list of every barcode and a list of every gene that our transcripts mapped to. My objective is to benchmark the operation of each tool and to compare the biological results obtained. To accompllish the latter, the chosen datasets are taken from the [Tabula Muris project](https://tabula-muris.ds.czbiohub.org/), a compendium of single cell transcriptome data from the model organism *Mus musculus* (The Tabula Muris Consortium, 2018). The count matrices obtained are then analyzed through an R notebook that seeks to emulate the process followed in the Tabula Muris study so that my results can be compared against those of the study.

## Progress as of 11/12/2022

I have completed all the pipelines, and I am running them with the trachea fastqs, available at GSM3040917. With the result obtained I am ironing out the details in the downstream notebook. The metrics with which to compare the results remain to be decided, and I am currently looking into a method to identify the cell types of the t-SNE clusters. I am also rewriting the analysis process provided by Tabula Muris, as it was written using an outdated version of the seurat package.

## Remaining tasks

* Bug fix and clean pipelines.

* Add scripts to generate accessory files to the git repository.

* Add description on how to obtain Tabula Muris raw data (relatively easy, using a specific 10X tool, but for completion and reproducibility).

* Decide on which metrics to compare the benchmark results (currently the pipelines produce a file from the bash time utility).

* Decide on which metrics to compare the downstream results.

* Produce code to present the comparisons.

* Analyze data from other tissues.

* Clean the notebook and improve the quality of the report.

* Extend this readme.

# Workflow and usage of the repository

## Preprocessing

The pipelines inside the pipeline folder take as input the paths of the "module.sh" file, of the raw data directory and of the mapping index, specific to each tool. Aditionally, a transcripts to genes file and/or tech-specific whitelist from 10X can be required, as well as other parameters. These inputs should be provided by directly editing the scripts. When the script is ready to run, do so in a terminal inside the directory where you want the output.

## Downstream

After the preprocessing is done, input the path of the three output files into the corresponding chunk in the notebook. In another chunk, the parameters used for the analysis can be specified. The defaults are the parameters used for the 8_15 trachea dataset from Tabula Muris. After this, the notebook can be knited to produce a report or run chunk by chunk.

# Bibliography

The Tabula Muris Consortium., Overall coordination., Logistical coordination. et al. Single-cell transcriptomics of 20 mouse organs creates a Tabula Muris. Nature 562, 367–372 (2018). https://doi.org/10.1038/s41586-018-0590-4
