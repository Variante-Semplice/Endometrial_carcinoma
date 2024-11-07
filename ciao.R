# ROMANEL PROJECT

#------------------------------
# Point 1: Load the RData file
#------------------------------
# per ketty (voi mettete il vostro path)
load("C:/Users/bello/Desktop/Endometrial_carcinoma/Uterine_corpus_endometrial_carcinoma.RData")
# raw_counts_df = contains the raw RNA-seq counts; 
# c_anno_df = contains sample names and conditions (case or control); 
# r_ anno_df = contains the ENSEMBL genes ids, the length of the genes and the genes symbols. 

#-----------------------------------------------
# Point 2: extracting only protein coding genes
#-----------------------------------------------
# we need to find which genes are coding, add a new column?
library(biomaRt)
# human genes dataset
ensembl <- useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl")
attributes <- listAttributes(ensembl)
# Find function for all genes 
g_f <- getBM(attributes=c("ensembl_gene_id", "gene_biotype"),mart = ensembl,
             filters = ("ensembl_gene_id"), values = list(r_anno_df$ensembl_gene_id))
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# mancano degli id ho solo 62034 elementi dei 62872
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# provo a cercare con gene name?
# ignorare sotto
g_f <- getBM(attributes=c("ensembl_gene_name", "gene_biotype"),mart = ensembl,
             filters = ("ensembl_gene_id"), values = list(r_anno_df$ensembl_gene_id))


g_coding <- g_f[which(g_f$gene_biotype == "protein_coding"),]

r_anno_df$use <- c(rep(NaN,length(r_anno_df$ensembl_gene_id)))

PC_r_anno <- r_anno[which(r_anno_df$ensembl_gene_id)]
# PC_raw_counts_df <-