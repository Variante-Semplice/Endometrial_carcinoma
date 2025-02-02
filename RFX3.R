# ROMANEL PROJECT

#-------------------------------
# Point 1: Load the RData file
#-------------------------------
load("Uterine_corpus_endometrial_carcinoma.RData")

#-----------------------------------------------
# Point 2: Extracting only protein coding genes
#-----------------------------------------------
# Use biomaRt package to retrieve the needed information; 
# Next tasks should use the new data-frames you have created. 

library(biomaRt)
# connecting to the human genes dataset
ensembl <- useEnsembl(biomart="ensembl",dataset="hsapiens_gene_ensembl", mirror = "asia")
# have a look at the possible attributes
attributes <- listAttributes(ensembl)
# Find function for all genes 
length(unique(r_anno_df$ensembl_gene_id)) 
# we have unique IDs in our dataset, good
g_f <- getBM(attributes = c("ensembl_gene_id", "gene_biotype", "external_gene_name"),
             mart = ensembl, filters = ("ensembl_gene_id"), 
             values = (r_anno_df$ensembl_gene_id))
# g_f contains ID, biotype and gene name of the genes in r_anno_df

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# mancano degli id ho solo 62034 elementi dei 62872
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# we try to find the missing genes with external name
geni_conId <-  r_anno_df$ensembl_gene_id %in%  g_f$ensembl_gene_id
# TRUE/FALSE indexes of genes that we found
r_geni_senzaId <- r_anno_df[!(geni_conId),]
# table of not found genes
# some do have the external gene name
# we try to find them in ensembl via  gene name
g_f1 <- getBM(attributes = c("ensembl_gene_id", "gene_biotype", "external_gene_name"),mart = ensembl,
              filters = ("external_gene_name"), values = (r_geni_senzaId$external_gene_name))
# we find 26 genes of the 838 missing but no one of this is gene coding
# so we continue with the ones that we have (62034)

g_coding <- g_f[which(g_f$gene_biotype == "protein_coding"),]
PC_r_anno <- r_anno_df[which(r_anno_df$ensembl_gene_id %in% g_coding$ensembl_gene_id),]
PC_raw_counts <- raw_counts_df[which(row.names(raw_counts_df) %in% g_coding$ensembl_gene_id),]
# we filter all the dataframes to keep only the coding genes

#--------------------------------------------
# Point 3: Differential Expression Analysis
#--------------------------------------------
# Perform a differential expression analysis using edgeR package and select up- 
# and down-regulated genes using an adjusted p-value cutoff of 0.01, a log fold 
# change ratio >1.5 for up-regulated genes and < (-1.5) for down-regulated genes 
# and a log CPM >1. Relax the thresholds if no or few results are available.  
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
# filter_vec is a vector with the genes and the number of replicates? for each
#Subset the raw counts matrix to include only genes that meet the filtering criteria

# check the dimension of the filtered matrix 
dim(filter_counts_df)
# we keep 25557 genes of the 62872

# apply the filter on gene annotation
filter_anno_df <- r_anno_df[rownames(filter_counts_df),] 
#Subset the annotation data to match the filtered raw counts data
dim(filter_anno_df)

#barplot of reads for each sample
# one column with the sample and one with total reads of that sample
size_df <- data.frame("sample"=colnames(filter_counts_df), 
                      "read_millions"=colSums(filter_counts_df)/1000000) 

ggplot(data=size_df,aes(sample,read_millions)) +
  geom_bar(stat="identity",fill="indianred2",colour="indianred2",width=0.7,alpha=0.7)+
  coord_flip()+
  theme_bw()

# boxplot of gene counts for each sample
long_counts_df <- gather(as.data.frame(filter_counts_df), key = "sample", value = "read_number")
# we put in a column the all the cells of the matrix and in the other the # of reads

ggplot(data=long_counts_df,aes(sample,read_number+1)) + 
  geom_boxplot(colour="cyan4",fill="cyan4",alpha=0.7) +
  theme_bw() +
  scale_y_log10()

##Prepare Data for Differential Expression Analysis

# create a DGRList object
edge_l <- DGEList(counts=filter_counts_df, group=c_anno_df$condition, 
                  samples=c_anno_df, genes=filter_anno_df) 
edge_l
# we put togheter reads, gene and sample infos

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
library(fgsea)
library(org.Hs.eg.db)
library(biomaRt)
library(clusterProfiler)
library(enrichplot)
library(ggnewscale)
library(DOSE)
library(pathview)
library(tidyverse)

### Load results of DEG (differentially expressed genes) analysis
DEGs <- read.table("DEGs.txt",header=T,sep="\t",as.is=T)
table(DEGs$class) #checks the count of genes by class (upregulated, downregulated or the same)

### Use biomaRt to map Gene symbols, Entrez IDs and Ensembl gen IDs
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "asia")
convert <- getBM(attributes=c("ensembl_gene_id","entrezgene_id","external_gene_name"),
                 filters=c("ensembl_gene_id"), 
                 values=DEGs$ensembl_gene_id,
                 mart = ensembl)
#Connects to the Ensembl database to retrieve entrezgene_id and 
#external_gene_name for each gene in DEGs.
#convert will store these mappings.
#merge() aligns DEG data with the convert mappings by ensembl_gene_id
DEGs <- merge(DEGs,convert)
DEGs <- DEGs[which(!is.na(DEGs$entrezgene_id)),]
DEGs <- DEGs[-which(duplicated(DEGs$entrezgene_id)),]
up_DEGs <- merge(up_DEGs,convert)
up_DEGs <- up_DEGs[which(!is.na(up_DEGs$entrezgene_id)),]
up_DEGs <- up_DEGs[-which(duplicated(up_DEGs$entrezgene_id)),]
down_DEGs <- merge(down_DEGs,convert)
down_DEGs <- down_DEGs[which(!is.na(down_DEGs$entrezgene_id)),]
down_DEGs <- down_DEGs[-which(duplicated(down_DEGs$entrezgene_id)),]

### GSEA analysis

## Create vector with ranks
ranks <- DEGs$logFC
ranks <- sort(ranks,decreasing = T)
names(ranks) <- DEGs$entrezgene_id
head(ranks)
barplot(sort(ranks, decreasing = T),las=3)
#ranks is created with log fold change (log2_FC) values from DEGs, 
#sorted in descending order.
#names(ranks) maps the ranks to entrezgene_id

### GO (Gene Ontology) analysis

## Perform Gene Ontology enrichment analysis (Biological Process)
#GO for BP on upregulated genes
ego_BP_up <- enrichGO(gene = up_DEGs$external_gene_name,
                      OrgDb = org.Hs.eg.db,
                      keyType = 'SYMBOL',
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)
#enrichGO() performs GO enrichment for Biological Process (ont = "BP") using upregulated genes.

#GO for BP on downregulated genes
ego_BP_down <- enrichGO(gene = down_DEGs$external_gene_name,
                        OrgDb = org.Hs.eg.db,
                        keyType = 'SYMBOL',
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05)
#show the top enriched GO terms and associated genes.

## Visualize the top 10 enriched terms with a barplot 
barplot(ego_BP_up,showCategory=10, main= "Top 10 enriched upregulated genes")
barplot(ego_BP_down,showCategory=10, main= "Top 10 enriched downregulated genes")

## Perform Gene Ontology (GO) enrichment analysis (Molecular Function)
#GO for MF on upregulated genes
ego_MF_up <- enrichGO(gene = up_DEGs$external_gene_name,
                      OrgDb = org.Hs.eg.db,
                      keyType = 'SYMBOL',
                      ont = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)
#GO for MF on downregulated genes
ego_MF_down <- enrichGO(gene = down_DEGs$external_gene_name,
                        OrgDb = org.Hs.eg.db,
                        keyType = 'SYMBOL',
                        ont = "MF",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05)

barplot(ego_MF_up,showCategory=10, main= "Top 10 enriched upregulated genes")
barplot(ego_MF_down,showCategory=10, main= "Top 10 enriched downregulated genes")

##WP analysis
#WP analysis, or WikiPathways (WP) analysis, is a type of pathway enrichment 
#analysis that uses pathways from the WikiPathways database. 
#WikiPathways is an open-access, community-driven database where scientists 
#can collaboratively curate and maintain biological pathways across different 
#species. It provides a wide range of pathways, including metabolic, signaling, 
#regulatory, and disease pathways, much like KEGG or Reactome but with the 
#added advantage of being open and editable.

#upregulated genes
eWP_up <- enrichWP(gene = up_DEGs$entrezgene_id,
                   organism = 'Homo sapiens',
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1)

head(eWP_up, n=10)

#downregulated genes
eWP_down <- enrichWP(gene = down_DEGs$entrezgene_id,
                     organism = 'Homo sapiens',
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.1)

head(eWP_up, n=10)


#-------------------------------------
# Point 5: Visualize one pathway
#-------------------------------------
## Perform KEGG enrichment analysis
eWP_KEGGS <- enrichKEGG(gene = up_DEGs$entrezgene_id,
                        organism = 'human',
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.1)
head(eWP_KEGGS, n=20)
#                   category              subcategory       ID
# hsa04110 Cellular Processes    Cell growth and death hsa04110
knitr::kable(head(eWP_KEGGS[, 1:6], n=10))

logFC <- up_DEGs$logFC
names(logFC) <- up_DEGs$entrezgene_id 
pathview(gene.data = logFC, 
         pathway.id = "hsa04110", 
         species = "human")

#------------------------------------
# Point 6: Identifying enriched TFs
#------------------------------------
library(biomaRt)
library(MotifDb)  # Access curated TF motif databases.
library(seqLogo)   # For visualizing DNA motifs using sequence logos.
library(PWMEnrich)  # Main package for performing motif enrichment analysis.
library(PWMEnrich.Hsapiens.background)

# finding promoter sequences
# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene") 
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg38") 

library(TxDb.Hsapiens.UCSC.hg38.knownGene)  # Transcript database for human genes (hg38 genome build).
library(BSgenome.Hsapiens.UCSC.hg38)       # Human genome sequence (hg38).

genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
ensembl <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "asia")
df <- getBM(attributes = c("external_gene_name",'entrezgene_id'),
            values=names(genes),filters ='entrezgene_id', mart = ensembl)
matched_indices <- match(genes$gene_id, df$entrezgene_id)
names(genes) <- ifelse(!is.na(matched_indices), df$external_gene_name[matched_indices], NA)

genes_name <- up_DEGs$external_gene_name # List of up-regulated genes
intersect_names <- intersect(names(genes), genes_name) #Filters genes to retain only those in both the up-regulated list and the UCSC database.

x <- promoters(genes,upstream = 500,downstream = 0)[intersect_names[intersect_names!= ""]] #Extracts promoter regions 500 bp upstream of the transcription start site (TSS).
seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38,x) 

# Calculating motif enrichment scores
library(MotifDb)   #Provides access to a curated collection of DNA sequence motifs.
library(seqLogo)   #Used for creating sequence logo plots.
library(PWMEnrich) #Offers functions for motif enrichment analysis.
library(PWMEnrich.Hsapiens.background)
data(PWMLogn.hg19.MotifDb.Hsap)
res = motifEnrichment(seq,PWMLogn.hg19.MotifDb.Hsap,score = "affinity")
report = sequenceReport(res, 1)
report
plot(report[report$p.value < 0.01], fontsize=7, id.fontsize=6)
plot(report[1:5])

#----------------------------------
# Point 7: PWM from MotifDB
#----------------------------------
# Select one among the top enriched TFs, compute the empirical distributions of 
# scores for all PWMs that you find in MotifDB for the selected TF and determine 
# for all of them the distribution (log2) threshold cutoff at 99.75% 
# (relax the threshold if needed)
# ho usato tp63 perchè sembrava un bel nome ma potremmo sceglierne altri
library(MotifDb) #Provides a curated collection of transcription factor binding motifs.
library(seqLogo)  #Used for visualization of sequence motifs as sequence logos.
library(PWMEnrich)  #Provides functionality for PWM-based motif enrichment analysis.
library(PWMEnrich.Hsapiens.background) #Supplies human-specific background models for enrichment.

mdb.human.selected = subset(MotifDb, organism=='Hsapiens' & geneSymbol=="RFX3")
if (length(mdb.human.selected) == 0) stop("No motifs found for the selected TF.")
motifs = as.list(mdb.human.selected)
length((mdb.human.selected))
pwms <- lapply(motifs, toPWM)
if (length(pwms) == 0) stop("No PWMs were generated for the selected TF.")
ecdf = motifEcdf(pwms,organism = "hg19",quick=TRUE)
# Initialize an empty list to store threshold values
thresholds <- numeric()
# Loop through each PWM's eCDF and calculate the threshold
for (ecdf_obj in ecdf) {
  # Calculate the 99.75 percentile threshold from the eCDF
  threshold <- log2(quantile(ecdf_obj, 0.9975))
  thresholds <- c(thresholds, threshold)  # Append to thresholds
}

# Display thresholds
threshold_data <- data.frame(TF = "RFX3", PWM = names(pwms), Threshold = thresholds)
print(threshold_data)

#funzionano ma fa caga (tra l'altro non ho capito perchè lo facciamo?)
scores = motifScores(seq,pwms,raw.score=TRUE)
plotMotifScores(scores,sel.motifs=names(mdb.human.selected),cols=c("red","green","blue"),cutoff=0.9975) 


#----------------------------------
# Point 8: Pattern Matching
#----------------------------------
# Identify which up-regulated genes have a region in their promoter 
# (defined as previously) with binding scores above the computed thresholds for 
# any of the previously selected PWMs. Use pattern matching as done during the course;

#non so se è giusto
scores = motifScores(seq, pwms, raw.scores = FALSE, verbose = FALSE, cutoff = thresholds)
highscore_seq <- which(apply(scores,1,sum)>0)
genes_id <- up_DEGs[highscore_seq,]

dim(genes_id)
head(genes_id)

#------------------------------------
# Point 9: PPI analysis
#------------------------------------
#Use STRING database to find PPI interactions among differentially expressed 
#genes and export the network in TSV format. 
library(igraph)
library(tidyr)
library(biomaRt)

write.table(unique(up_DEGs$external_gene_name),sep = '\t', file = 'PPI_up_DEGs.txt',row.names = F, 
            col.names = F, quote = T)

links <- read.delim("string_interactions_short.tsv") # TSV data downloaded from STRING webpage 
#(dovete scaricarlo da STRING perchè è troppo pesante per git)
head(links)
UP__DEGs <- read.table("PPI_up_DEGs.txt", sep = "\t", header = TRUE)


#--------------------
# Point 10: Network
#--------------------
## Create nodes annotations using biomaRt
library(biomaRt)
library(igraph)
library(dplyr)

# Step 1: Use biomaRt to retrieve node annotations
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = 'asia')

nodes <- getBM(attributes = c("external_gene_name", "ensembl_gene_id", "description", 
                              "gene_biotype", "start_position", "end_position", 
                              "chromosome_name", "strand"),
               filters = "external_gene_name", 
               values = unique(c(links$X.node1, links$node2)),  # Combine all nodes in the network
               mart = ensembl)


nodes <- getBM(attributes=c("external_gene_name","ensembl_gene_id","description","gene_biotype","start_position","end_position","chromosome_name","strand"),
               filters=c("external_gene_name"), 
               values=UP__DEGs [,1],
               mart = ensembl)




# Keep only relevant columns and remove duplicates
nodes <- unique(nodes[, c("external_gene_name", "description", "gene_biotype", 
                          "start_position", "end_position", "chromosome_name", "strand")])

# Step 2: Filter out nodes that are not in the links
valid_nodes <- unique(c(links$X.node1, links$node2))  # All nodes in the links
## Remove duplicate entries in the nodes data frame
# Ensure external_gene_name is unique
nodes <- nodes %>%
  distinct(external_gene_name, .keep_all = TRUE)

# Build the graph after deduplication
net <- graph_from_data_frame(
  d = links[, c("X.node1", "node2")],  # Edges
  vertices = nodes,                   # Filtered unique nodes
  directed = FALSE
)

# Check the network object
class(net)  # Should be "igraph"
net


## Check attributes
edge_attr(net)
vertex_attr(net)
graph_attr(net)

## Plot the PPI network
plot(net, 
     edge.width=2,
     vertex.color="orange",
     vertex.size=10,
     vertex.frame.color="darkgray",
     vertex.label.color="black", 
     vertex.label.cex=0.7,
     edge.curved=0.1) 




# Find the node with the highest degree (most connections)
most_connected_node_23 <- V(net)[which.max(degree(net))]

# Create subgraph of this node and its immediate neighbors
neighborhood_net_23 <- make_ego_graph(net, order = 1, nodes = most_connected_node_23)[[1]]
plot(neighborhood_net_23)