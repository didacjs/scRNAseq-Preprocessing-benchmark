library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl") # get mouse biomart
gene_description <-getBM(attributes = c('external_gene_name', 'description'), # get table with ensembl id and description
                         mart = ensembl)
names(gene_description)[1] <- "gene" # change name to match the ensembl column in the top table
write.csv(gene_description, "~/Work/Code/scRNAseq-Preprocessing-benchmark/Refs/extid2description.csv", row.names=FALSE) # write to file
ensemblid2description<-read_csv("~/Work/Code/scRNAseq-Preprocessing-benchmark/Refs/extid2description.csv") # check file
