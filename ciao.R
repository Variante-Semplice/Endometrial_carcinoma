# ROMANEL PROJECT

#-------------------------------
# Point 1: Load the RData file
#-------------------------------
# per ketty (voi mettete il vostro path)
load("Uterine_corpus_endometrial_carcinoma.RData")
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
length(unique(r_anno_df$ensembl_gene_id)) # unique IDs
g_f <- getBM(attributes = c("ensembl_gene_id", "gene_biotype", "external_gene_name"),mart = ensembl,
             filters = ("ensembl_gene_id"), values = (r_anno_df$ensembl_gene_id))

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# mancano degli id ho solo 62034 elementi dei 62872
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
geni_conId <-  r_anno_df$ensembl_gene_id %in%  g_f$ensembl_gene_id
r_geni_senzaId <- r_anno_df[!(geni_conId),]
# alcuni hanno external gene name
# li cerchiamo mediante gene name
g_f1 <- getBM(attributes = c("ensembl_gene_id", "gene_biotype", "external_gene_name"),mart = ensembl,
              filters = ("external_gene_name"), values = (r_geni_senzaId$external_gene_name))
# we find 26 genes of the 838 missing but no one of this is gene coding
# so we continue with the ones that we have (62034)

g_coding <- g_f[which(g_f$gene_biotype == "protein_coding"),]
PC_r_anno <- r_anno_df[which(r_anno_df$ensembl_gene_id %in% g_coding$ensembl_gene_id),]
PC_raw_counts <- raw_counts_df[which(row.names(raw_counts_df) %in% g_coding$ensembl_gene_id),]

#--------------------------------------------
# Point 3: differential Expression Analysis
#--------------------------------------------
library("GenomicFeatures")
library("ggplot2")
library("stringr")
library("tidyverse")
library("edgeR")
library("pheatmap")

## Filter raw counts data retaining only genes with a raw count >20 in at least 
## 5 Cases or 5 Control samples

count_thr <- 20 #minimum count for a gene to be considered present
repl_thr <- 5 #ensures that we retain genes with counts above 20 in at least 5 samples in either group.

filter_vec <- apply(raw_counts_df,1,
                    function(y) max(by(y, c_anno_df$condition, function(x) sum(x>=count_thr))))
# this vector checks if the gene meets the threshold in either condition 
# by grouping counts by condition and checking if there are at least 
# repl_thr samples with counts above count_thr
# we have 30 case and 30 control in c_anno_df, there are some genes which have
# more than 20 reads for all the 30 case/control.

#statistics for the filtering
table(filter_vec)

filter_counts_df <- raw_counts_df[filter_vec>=repl_thr,] 
#Subset the raw counts matrix to include only genes that meet the filtering criteria

# check the dimension of the filtered matrix 
dim(filter_counts_df)
# we keep 25557 genes

# apply the filter on gene annotation
filter_anno_df <- r_anno_df[rownames(filter_counts_df),] 
#Subset the annotation data to match the filtered raw counts data
dim(filter_anno_df)

#barplot of reads for each sample
size_df <- data.frame("sample"=colnames(filter_counts_df), 
                      "read_millions"=colSums(filter_counts_df)/1000000) 

ggplot(data=size_df,aes(sample,read_millions)) +
  geom_bar(stat="identity",fill="indianred2",colour="indianred2",width=0.7,alpha=0.7)+
  coord_flip()+
  theme_bw()

#boxplot of gene counts for each sample
long_counts_df <- gather(as.data.frame(filter_counts_df), key = "sample", value = "read_number")

ggplot(data=long_counts_df,aes(sample,read_number+1)) + 
  geom_boxplot(colour="cyan4",fill="cyan4",alpha=0.7) +
  theme_bw() +
  scale_y_log10()

##Prepare Data for Differential Expression Analysis

# create a DGRList object
edge_l <- DGEList(counts=filter_counts_df, group=c_anno_df$condition, samples=c_anno_df, genes=filter_anno_df) 
edge_l

# normalization
edge_n <- calcNormFactors(edge_l, method="TMM")
edge_n

# create a cpm table (normalized expression values)
cpm_table <- as.data.frame(round(cpm(edge_n),2))
head(cpm_table)

#boxplot after normalization
long_cpm_df <- gather(cpm_table, key = "sample", value = "CPM")

ggplot(data=long_cpm_df,aes(sample,CPM+1)) +
  geom_boxplot(colour="orchid3",fill="orchid3",alpha=0.7)+
  theme_bw()+
  scale_y_log10()

group <- edge_n$samples$group 
# define the experimental design matrix
design <- model.matrix(~0+group, data=edge_n$samples)
colnames(design) <- levels(edge_n$samples$group)
rownames(design) <- edge_n$samples$sample
design

# calculate dispersion and fit with edgeR
edge_d <- estimateDisp(edge_n, design)
edge_f <- glmQLFit(edge_d, design) 

# definition of the contrast (conditions to be compared)
contro <- makeContrasts("case - control", levels = design)

# fit the model with generalized linear models
edge_t <- glmQLFTest(edge_f,contrast=contro)
DEGs <- as.data.frame(topTags(edge_t,n=20000,p.value = 0.01, sort.by = "logFC"))
DEGs$class <- "=" # Initialize a 'class' column to categorize genes
DEGs$class[which(DEGs$logCPM > 1 & DEGs$logFC > 1.5 & DEGs$FDR < 0.01)] <- "+"  # Up-regulated
DEGs$class[which(DEGs$logCPM > 1 & DEGs$logFC < (-1.5) & DEGs$FDR < 0.01)] <- "-"  # Down-regulated

# Sort the DEGs by log fold-change in descending order
DEGs <- DEGs[order(DEGs$logFC,decreasing = T),]

head(DEGs)
table(DEGs$class)

## Volcano plot
# Add a column for -log10(FDR) for plotting purposes
DEGs$neg_log10_FDR <- -log10(DEGs$FDR)

# Create a volcano plot 
volcano_plot <- ggplot(DEGs, aes(x = logFC, y = neg_log10_FDR, color = class)) +
  geom_point(alpha = 0.6, size = 1.5) +  # Points with transparency and size
  scale_color_manual(values = c("=" = "grey", "+" = "tomato3", "-" = "forestgreen")) + 
  labs(
    title = "Volcano Plot of Differential Expression",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-Value (FDR)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14)
  ) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "black") +  # FDR threshold line
  geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed", color = "black")    # Fold-change thresholds

print(volcano_plot)

##Heatmap
cols <- cols <- c(ifelse(c_anno_df$condition == "case", "sienna2", "goldenrod1")) 
pal <- c("forestgreen","white","tomato3") 
pal <- colorRampPalette(pal)(100)

heatmap(as.matrix(cpm_table[which(rownames(cpm_table)%in%DEGs$ensembl_gene_id[which(DEGs$class!="=")]),]),
        ColSideColors = cols,cexCol = 0.5,margins = c(4,4),col=pal,cexRow = 0.2)

## Export differentially expressed genes in a text file
up_DEGs <- DEGs[which(DEGs$class=="+"),]
down_DEGs <- DEGs[which(DEGs$class=="-"),]

write.table(up_DEGs,file="up_DEGs.txt",row.names=F,col.names=T,sep="\t",quote=F)
write.table(down_DEGs,file="down_DEGs.txt",row.names=F,col.names=T,sep="\t",quote=F)
write.table(DEGs,file="DEGs.txt",row.names=F,col.names=T,sep="\t",quote=F)

#-------------------------------------
# Point 4: GSEA using clusterProfiler
#-------------------------------------
##Gene Set Enrichment Analysis (GSEA) is a statistical technique used to 
##determine whether a predefined set of genes shows statistically significant, 
##coordinated differences in expression under two different conditions, 
##such as "treatment" vs. "control" or "healthy" vs. "disease." 
##Instead of evaluating individual genes separately, GSEA assesses entire groups 
##of genes that are functionally related, such as those involved in a particular 
##pathway or biological process.