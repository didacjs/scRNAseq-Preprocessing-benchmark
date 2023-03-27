# Benchmark and downstream comprative of four scRNA-Seq preprocessing pipelines (W.I.P.)

## General description

Single cell transcriptomic techniques (scRNA-Seq) examine the gene expression level in individual cells, simultaneously measuring concentration of messenger RNA transcripts (mRNA) belonging to hundreds to thousands of genes. The analysis of transcriptomics at this  resolution makes possible the decryption of heterogeneous cell populations, the study of cellular life cycles and the modeling of transcriptional dynamics.

The raw output of these techniques are sequencing files (FASTQ) that contain millions of sequences (reads) of transcripts. Each read contains a nucleotide sequence encoding a barcode and a Unique Molecular Identifier (UMI). The barcode associates the read with a cell and the UMI with an RNA molecule. These files cannot be analyzed directly. As a first step, it is necessary to preprocess the raw data to produce a counts matrix, which contains the number of  transcripts for each gene in each cell. With the counts matrix, we can then study the transcriptomic profile of our cell sample. 

The preprocessing step can be computationally intensive depending on the size of the sample, and so a variety of specialized tools to carry and optimize the task have been developed. Despite that, there is still no procedure considered as a standard. Consequently, it is interesting to compare the different available tools to find the one that best suits our needs, computational resources and specifications of the experiment.

## Approach

I have designed four scRNA-seq preprocessing pipelines in bash, each starting from one different tool, STAR/UMI-Tools, STARsolo, Kallisto|Bustools and Salmon Alevin. These pipelines take FASTQ files as input and generate a counts matrix, along with a list of every barcode and a list of every gene that our transcripts mapped to. My objective is to benchmark the operation of each tool and to compare the biological results obtained. To accompllish the latter, the chosen datasets are taken from the [Tabula Muris project](https://tabula-muris.ds.czbiohub.org/), a compendium of single cell transcriptome data from the model organism *Mus musculus* (The Tabula Muris Consortium, 2018). The count matrices obtained are then analyzed to compare each preprocessing tool.

I have chosen UMI-Tools/STAR because it was one of the first tools available for preprocessing single cell and will serve as a reference for other cutting-edge tools (Smith et al., 2017). STARsolo is an optimized workflow integration of the STAR aligner (Kaminow et al., 2021). On the other hand, Kallisto|Bustools and SalmonAlevin use pseudomapping to reduce execution time and the resources associated with the alignment of the reads to the genome (Melsted et al., 2021; Srivastava et al.,2019).

With these tools I have written pipelines as Bash scripts. The inputs of the  pipelines, in addition to the dataset (FASTQ reads), are the genome of the species and its annotation file. As input I will use the raw data from the Tabula Muris project. The main output will be a matrix whose format may vary, but which must indicate the count of transcripts of each gene in each cell. On the one hand, the scripts save data on the resources used and the runtime that I can use to compare the pipelines. On the other hand, I can compare the number of cells and the number of transcripts per cell of the matrices obtained with those of the matrices from Tabula Muris.

Finally, the downstream analysis will be carried out in a notebook in R, using the Seurat library (Hao et al., 2021). The downstream analysis can go much deeper, but I will focus on quality control, dimensional reduction, clustering and extraction of differentially expressed genes and lastly, annotation. Obtaining the most differentially expressed genes from each cluster and annotating them, my intention is to deduce what kind of functions are carried out by the type of cells in each cluster. I can compare the results with those obtained by doing the same analysis on the matrices provided by Tabula Muris. With this last comparison, I will be able to evaluate the biological results of the preprocessing tools.

## Updates

### Progress as of 27/03/2023

I implemented a way to compare the t-SNE clusters obtained from our pipelines to the same clusters obtained from already-preprocessed data from Tabula Muris. For this feature another .rmd is added as well as scripts to obtain supplementary files from BioMart.

### Progress as of 11/12/2022

I have completed all the pipelines, and I am running them with the trachea fastqs, available at GSM3040917. With the result obtained I am ironing out the details in the downstream notebook. The metrics with which to compare the results remain to be decided, and I am currently looking into a method to identify the cell types of the t-SNE clusters. I am also rewriting the analysis process provided by Tabula Muris, as it was written using an outdated version of the seurat package.

### Remaining tasks

* Bug fix and clean pipelines.

* Add scripts to generate accessory files to the git repository.

* Add description on how to obtain Tabula Muris raw data (relatively easy, using a specific 10X tool, but for completion and reproducibility).

* Decide on which metrics to compare the benchmark results (currently the pipelines produce a file from the bash time utility).

* Decide on which metrics to compare the downstream results.

* Analyze data from other tissues.

* Clean the notebook and improve the quality of the report.

* Extend this readme.

# Workflow and usage of the repository

## Preprocessing

The pipelines inside the pipeline folder take as input the paths of the "module.sh" file, of the raw data directory and of the mapping index, specific to each tool. Aditionally, a transcripts to genes file and/or tech-specific whitelist from 10X can be required, as well as other parameters. These inputs should be provided by directly editing the scripts. When the script is ready to run, do so in a terminal inside the directory where you want the output.

## Downstream

After the preprocessing is done, input the path of the three output files into the corresponding chunk in the notebook. In another chunk, the parameters used for the analysis can be specified. The defaults are the parameters used for the 8_15 trachea dataset from Tabula Muris. After this, the notebook can be knited to produce a report or run chunk by chunk.

# Bibliography

Hao, Y., Hao, S., Andersen-Nissen, E., Mauck, W. M., Zheng, S., Butler, A., Lee, M. J., Wilk, A. J., Darby, C.,Zager, M., Hoffman, P., Stoeckius, M., Papalexi, E., Mimitou, E. P., Jain, J., Srivastava, A., Stuart, T.,Fleming, L. M., Yeung, B., ... Satija, R. (2021). Integrated analysis of multimodal single-cell data. Cell,
184(13), 3573-3587.e29. https://doi.org/10.1016/J.CELL.2021.04.048/ATTACHMENT/1E5EB5C1-59EE-4B2B-8BFA-14B48A54FF8F/MMC3.XLSX
Kaminow, B., Yunusov, D., & Dobin, A. (2021). STARsolo: accurate, fast and versatile mapping/quantification of single-cell and single-nucleus RNA-seq data. BioRxiv,
2021.05.05.442755. https://doi.org/10.1101/2021.05.05.442755 

Melsted, P., Booeshaghi, A. S., Liu, L., Gao, F., Lu, L., Min, K. H. (Joseph), da Veiga Beltrame, E. Hjörleifsson, K. E., Gehring, J., & Pachter, L. (2021). Modular, Efficient and constant-memory single-cell RNA-seq preprocessing. Nature Biotechnology 2021 39:7, 39(7), 813–818.https://doi.org/10.1038/s41587-021-00870-2

Schaum, N., Karkanias, J., Neff, N. F., May, A. P., Quake, S. R., Wyss-Coray, T., Darmanis, S., Batson, J., Botvinnik, O., Chen, M. B., Chen, S., Green, F., Jones, R. C., Maynard, A., Penland, L., Pisco, A. O., Sit, R. V., Stanley, G. M., Webber, J. T., ... Weissman, I. L. (2018). Single-cell transcriptomics of 20
mouse organs creates a Tabula Muris. Nature 2018 562:7727, 562(7727), 367–372. https://doi.org/10.1038/s41586-018-0590-4

Smith, T., Heger, A., & Sudbery, I. (2017). UMI-tools: modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy. Genome Research, 27(3), 491–499. https://doi.org/10.1101/GR.209601.116

Srivastava, A., Malik, L., Smith, T., Sudbery, I., & Patro, R. (2019). Alevin efficiently estimates accurate gene abundances from dscRNA-seq data. Genome Biology, 20(1), 1–16. https://doi.org/10.1186/S13059-019-1670-Y/FIGURES/8
