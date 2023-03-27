library(biomaRt)
mart <- useEnsembl('ensembl', dataset = 'mmusculus_gene_ensembl')
rRNA<-biomaRt::getBM(values="rRNA", 
               filters="biotype", 
               attributes=c("ensembl_gene_id"), 
               mart = mart)
write_csv(rRNA, "~/Work/rRNA.csv", col_names = T)

