##--------------------------------------- Influenza biomarkers project overview
# The goal of this project is to computationally identify biomkarers that validate the efficacy of an in-vitro "organ-on-chip" model of influenza infection.
# In particular, we hope to identify biomarkers that capture the idiosyncracies of COPD -- a blanket term used to describe a variety of different lung diseases. 
# Obtaining these biomakers will aid in i) evaluating the model's ability to recapitulate disease physiology ii) identifying novel clinical biomarkers that 
# distinguish influenza from other respiratory diseases and iii) aid in the development of therapeutics tageting influenza. These biomarkers will be recovered 
# from RNA-sequencing data capturing gene expression dynamics in a healthy and infected (influenza) model.  This project is being conducted in collaboration with the 
# Ingber lab and will utilize their acquired RNA-seq data alongside publicly available datasets. 

# RNA-seq is of lung epithelial cells
# there are 5 COPD donors and 4 healthy donors
# for most conditions, we have 4hr time points of control vs. virus and 18hr time points of control, poly IC, and
# virus. Only exception is one healthy donor which only has three 18hr virus replicates and one healthy donor with 
# only four control conditions at 4 hrs

# Author: Miguel A. Alcantar 
# Last updated: 09/25/2019

##--------------------------------------- Loading libraries needed for analysis

# If starting from scratch, you will likely need to install the following packages using BiocManager. 
#if (!requireNamespace("BiocManager")) --- likely won't need to run this
#install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install("apeglm")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("DiffNet")
# BiocManager::install("hgu133plus2.db")
# BiocManager::install("limma")
# BiocManager::install("GEOquery")
# BiocManager::install("affy")
# BiocManager::install("ArrayExpress")
# BiocManager::install("arrayQualityMetrics")
# remotes::install_github("diogocamacho/rpegeos") # download rpegeos from Diogo's github
# devtools::install_github('kevinblighe/EnhancedVolcano')
# BiocManager::install("mygene")
# devtools::install_github("sinhrks/ggfortify")


# loading packages
# general packages for data organization and parsing files
library(tidyverse) # for organizing data structures
library(readr) # for data parsing / loading in data files
library(readxl) # for reading Excel files
library(EnhancedVolcano)
# library(gdtools)

# differential expression
library(DESeq2) # differential expression (RNA-seq)
library(limma, quietly = TRUE) # differential expression (microarray)

# gene database
library(org.Hs.eg.db) # gene database
library(DiffNet)
library(hgu133plus2.db, quietly = TRUE) # human genome

# pathway analysis
library(pathview)
library(gage)
library(gageData)
library(rpegeos) # created by Diogo (https://rdrr.io/github/diogocamacho/rpegeos/)
library(ggfortify)

data(kegg.sets.hs) # KEGG pathways
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
data(kegg.gs)


#library(DRUID) created by Diogo -- don't need yet (https://github.com/diogocamacho/druid/)
#library(grove)

##--------------------------------------- Loading scripts with functions used throughout analysis

source("util/load_transcriptomics_data_functions.R") # script for loading data
source("util/diff_expression_pathway_analysis_functions.R") # script for gene- and pathway-level analysis
source("util/preprocess_transcriptomics_data_functions.R") # scripts for preprocessing data

##--------------------------------------- Loading data and reorganizing to make easier to handle
# load metadata and add sample IDs that match count matrix
copd_nhb_RNAseq_original <- load_new_COPD_NHB_RNAseq_data()

copd_nhb_count_data_original <- copd_nhb_RNAseq_original$count_data
copd_nhb_genes_original <- copd_nhb_RNAseq_original$gene_ids
sample_data_copd_nhb <- copd_nhb_RNAseq_original$metadata # metadata

# final gene count and ID data matrices
copd_nhb_data_processed <- preprocess_chip_data(copd_nhb_count_data_original, copd_nhb_genes_original, org.Hs.eg.db, freq_thresh=0.50)
count_data_copd_nhb<- copd_nhb_data_processed$count_data # preprocessed count data 
gene_ids_copd_nhb <- copd_nhb_data_processed$gene_ids # genes remaining after preprocessing
genes_filtered <- copd_nhb_data_processed$genes_filtered
gkey_all <- AnnotationDbi::select(x = org.Hs.eg.db, 
                                  keys = gene_ids_copd_nhb, 
                                  columns = c("ENTREZID","GENENAME"), 
                                  keytype = "SYMBOL")  
# write.csv(sample_data_copd_nhb, file = "../../../data/RNA-seq_data_chip/COPD_NHB_data/sample_data_copd_nhb_preprocessed.csv")
# write.csv(sample_data_copd_nhb, file = "../../../data/RNA-seq_data_chip/COPD_NHB_data/sample_data_copd_nhb_preprocessed.csv")
write.csv(gkey_all, file = "../../../data/RNA-seq_data_chip/COPD_NHB_data/gkey_all_preprocessed.csv")

# write.csv(count_data_copd_nhb, file = "../../../data/RNA-seq_data_chip/COPD_NHB_data/count_data_copd_nhb_preprocessed.csv")
# write.table(gene_ids_copd_nhb, file= "../../../data/gene_lists/RNAseq_NHB_COPD/gene_ids_copd_nhb_preprocessed.txt", sep="\t", col.names = F, row.names = F, quote=F)

##--------------------------------------- Differential expression analysis with DESeq2

# differential expression ----
# use deseq2
dds_copd_nhb <- DESeqDataSetFromMatrix(countData = count_data_copd_nhb,
                                   colData = sample_data_copd_nhb,
                                   design = ~ condition_group_time)

dds_copd_nhb <- DESeq(object = dds_copd_nhb)
counts_copd_nhb_normalized <- counts(dds_copd_nhb, normalized=TRUE) # normalize counts using sizeFactor; useful for GSEA_GAGE


count_data_copd_nhb_for_pca <- as_tibble(scale(t()))%>% 
  dplyr::mutate(., condition_group_time = sample_data_copd_nhb$condition_group_time,
                condition_group = sample_data_copd_nhb$condition_group)

autoplot(prcomp(count_data_copd_nhb_for_pca[,1:(length(count_data_copd_nhb_for_pca)-2)], ), data = count_data_copd_nhb_for_pca, colour = 'condition_group')

# regular copd vs control
# compare basal differences between copd and nhb
contrast_copd_nhb_contr_4 <- c("condition_group_time", "COPD_control_4", "NHB_control_4")
contrast_copd_nhb_contr_18 <- c("condition_group_time", "COPD_control_18", "NHB_control_18")

# compare copd patients before and after viral infection
contrast_copd_v4 <- c("condition_group_time", "COPD_virus_4", "COPD_control_4")
contrast_copd_v18 <- c("condition_group_time", "COPD_virus_18", "COPD_control_18")

# compare nhb patients before and after viral infection
contrast_nhb_v4 <- c("condition_group_time", "NHB_virus_4", "NHB_control_4")
contrast_nhb_v18 <- c("condition_group_time", "NHB_virus_18", "NHB_control_18")

# compare copd and nhb patients with viral infections 
contrast_copd_nhb_v4 <- c("condition_group_time", "COPD_virus_4", "NHB_virus_4")
contrast_copd_nhb_v18 <- c("condition_group_time", "COPD_virus_18", "NHB_virus_18")

# compare copd and nhb patients with viral infections 
contrast_copd_poly_ic18 <- c("condition_group_time", "COPD_poly_ic_18", "COPD_control_18")
contrast_nhb_poly_icv8 <- c("condition_group_time", "NHB_poly_ic_18", "NHB_control_18")


res_copd_nhb_contr4 <- compare_group_DESeq2(dds_copd_nhb, contrast_copd_nhb_contr_4, gene_ids_copd_nhb, MHT = "fdr")
res_copd_nhb_contr18<- compare_group_DESeq2(dds_copd_nhb, contrast_copd_nhb_contr_18, gene_ids_copd_nhb, MHT = "fdr")

res_copd_v4<- compare_group_DESeq2(dds_copd_nhb, contrast_copd_v4, gene_ids_copd_nhb, MHT = "fdr")
res_copd_v18<- compare_group_DESeq2(dds_copd_nhb, contrast_copd_v18, gene_ids_copd_nhb, MHT = "fdr")

res_nhb_v4<- compare_group_DESeq2(dds_copd_nhb, contrast_nhb_v4, gene_ids_copd_nhb, MHT = "fdr")
res_nhb_v18<- compare_group_DESeq2(dds_copd_nhb, contrast_nhb_v18, gene_ids_copd_nhb, MHT = "fdr")

res_copd_nhb_v4<- compare_group_DESeq2(dds_copd_nhb, contrast_copd_nhb_v4, gene_ids_copd_nhb, MHT = "fdr")
res_copd_nhb_v18<- compare_group_DESeq2(dds_copd_nhb, contrast_copd_nhb_v18, gene_ids_copd_nhb, MHT = "fdr")

res_copd_poly_ic18<- compare_group_DESeq2(dds_copd_nhb, contrast_copd_poly_ic18, gene_ids_copd_nhb, MHT = "fdr")
res_nhb_poly_ic18<- compare_group_DESeq2(dds_copd_nhb, contrast_nhb_poly_icv8, gene_ids_copd_nhb, MHT = "fdr")


### DGCR5 barplots--------------------------
DGCR5_NHB_infected_4 <- res_nhb_v4$log2FoldChange[which(res_nhb_v4$gene == "DGCR5")]
DGCR5_NHB_infected_18 <- res_nhb_v18$log2FoldChange[which(res_nhb_v18$gene == "DGCR5")]

DGCR5_COPD_infected_4 <- res_copd_v4$log2FoldChange[which(res_copd_v4$gene == "DGCR5")]
DGCR5_COPD_infected_18 <- res_copd_v18$log2FoldChange[which(res_copd_v18$gene == "DGCR5")]

DGCR5_COPD_NHB_noninfected_4 <- res_copd_nhb_contr4$log2FoldChange[which(res_copd_nhb_contr4$gene == "DGCR5")]
DGCR5_COPD_NHB_noninfected_18 <- res_copd_nhb_contr18$log2FoldChange[which(res_copd_nhb_contr18$gene == "DGCR5")]

DGCR5_COPD_NHB_infected_4 <- res_copd_nhb_v4$log2FoldChange[which(res_copd_nhb_v4$gene == "DGCR5")]
DGCR5_COPD_NHB_infected_18 <- res_copd_nhb_v18$log2FoldChange[which(res_copd_nhb_v18$gene == "DGCR5")]

conditions = c("NHB infected (4hr)","NHB infected (18hr)",
               "COPD infected (4hr)", "COPD infected (18hr)", 
               "COPD_NHB non-infected (4hrs)","COPD_NHB non-infected (18hrs)",
               "COPD_NHB infected (4hrs)","COPD_NHB infected (18hrs)")

DBGR5_fold_change <- c(DGCR5_NHB_infected_4,DGCR5_NHB_infected_18,
                       DGCR5_COPD_infected_4,DGCR5_COPD_infected_4,
                       DGCR5_COPD_NHB_noninfected_4,DGCR5_COPD_NHB_noninfected_18,
                       DGCR5_COPD_NHB_infected_4, DGCR5_COPD_NHB_infected_18)

DGCR5_foldchange_df <- data.frame(condition=conditions,
                 fold_change=DBGR5_fold_change)
DGCR5_foldchange_df$condition <-factor(DGCR5_foldchange_df$condition
                                       , levels = unique(conditions))

p<-ggplot(data=DGCR5_foldchange_df, aes(x=condition, y=fold_change)) + geom_bar(stat="identity") + 
 theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                    axis.text.x = element_text(angle = 45, hjust = 1, size=14),
                    axis.text.y = element_text( size=14),
                    axis.title.x = element_text(size = 14),
                    axis.title.y = element_text(size = 14),
                    plot.title = element_text(size = (14), hjust=0.5)) + labs(x = "Condition", y = "log2 fold change", title= "DGCR5 lncRNA log2 fold changes")
p

ggsave(filename = "DGCR5_lncRNA_foldchange.pdf",
       plot = last_plot(),
       device = "pdf",
       height = 8,
       width = 10,
       scale = 1,
       dpi = 400
)

### DGCR5 barplots raw transcript counts--------------------------
DGCR5_df <- data.frame() #one column for transcript count; one column for 
conditions <- unique(sample_data_copd_nhb$condition_group_time)
counter <- 1
DGCR5_idx <- which(gene_ids_copd_nhb == "DGCR5")
for(condition in conditions){
  condition_idx <- which(sample_data_copd_nhb$condition_group_time == condition)
  for(idx in condition_idx){
    DGCR5_df[counter, "Condition"] <- condition
    DGCR5_df[counter, "transcript_count"] <- count_data_copd_nhb[DGCR5_idx,idx]
    counter <- counter+1
  }
}
DGCR5_df$Condition <- factor(DGCR5_df$Condition, levels = c("NHB_control_4","NHB_control_18",
                                                            "COPD_control_4","COPD_control_18",
                                                            "NHB_virus_4" ,"NHB_virus_18",
                                                            "COPD_virus_4","COPD_virus_18",
                                                            "NHB_poly_ic_18","COPD_poly_ic_18"
                                                            ),ordered = TRUE)

my_comparisons <-list(c("NHB_virus_4",  "NHB_virus_18"),c("COPD_virus_4",  "COPD_virus_18"),
                      c("NHB_control_4",  "NHB_control_18"),c("COPD_control_4",  "COPD_control_18"),
                      c("COPD_control_4",  "NHB_control_4"),
                      c("COPD_control_18",  "NHB_control_18"),
                      c("COPD_virus_4",  "NHB_virus_4"),
                      c("COPD_virus_18",  "NHB_virus_18"))
ggplot(DGCR5_df, aes_string(x = "Condition", y="transcript_count")) + geom_boxplot()+
  theme_classic()+stat_compare_means(comparisons=my_comparisons,
                                    label = "p.format",
                                    method = "t.test",
                                    paired=FALSE, method.args = list(var.equal = FALSE))+
  labs(x="Condition" ,y="transcript counts")+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                                  axis.text.x = element_text(angle = 45, hjust = 1, size=14),
                                                  axis.text.y = element_text( size=14),
                                                  axis.title.x = element_text(size = 14),
                                                  axis.title.y = element_text(size = 14),
                                                  plot.title = element_text(size = (14), hjust=0.5))

ggsave(filename = "../../../figs/lncRNAs/DGCR5_lncRNA_transcript_counts.pdf",
       plot = last_plot(),
       device = "pdf",
       height = 8,
       width = 10,
       scale = 1,
       dpi = 400
)

### DGCR5 correlation--------------------------
DGCR5_idx <- which(copd_nhb_genes_original == "DGCR5")

DGCR5_counts<- copd_nhb_count_data_original[DGCR5_idx,]
IFNB_idx <- which(copd_nhb_genes_original == "IFNB1")
IFNB_counts<- copd_nhb_count_data_original[IFNB_idx,]

cor(DGCR5_counts, IFNB_counts, method = "pearson")
DGCR5_IFNB_df <- data.frame(matrix(nrow = 162))
DGCR5_IFNB_df$DGCR5 <- as.matrix(DGCR5_counts)[1,]
DGCR5_IFNB_df$IFNB <- as.matrix(IFNB_counts)[1,]

ggscatter( DGCR5_IFNB_df, x ="DGCR5", y = "IFNB", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "DGCR5 transcript count", ylab = "IFNB transcript count", title = "IFNB vs. DGCR5 all samples")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
               axis.text.x = element_text( size=14),
               axis.text.y = element_text( size=14),
               axis.title.x = element_text(size = 14),
               axis.title.y = element_text(size = 14),
               plot.title = element_text(size = (14), hjust=0.5))

ggsave(filename = "../../../figs/lncRNAs/DGCR5_IFNB_correlate_allsamples.pdf",
       plot = last_plot(),
       device = "pdf",
       height = 8,
       width = 10,
       scale = 1,
       dpi = 400
)
COPD_virus_18_idx <- which(sample_data_copd_nhb$condition_group_time == "COPD_virus_18")
NHB_virus_18_idx <- which(sample_data_copd_nhb$condition_group_time == "NHB_virus_18")

DGCR5_IFNB_COPD_virus_18_df <- DGCR5_IFNB_df[COPD_virus_18_idx,]
ggscatter( DGCR5_IFNB_COPD_virus_18_df, x ="DGCR5", y = "IFNB", 
           add = "reg.line", conf.int = TRUE, 
           cor.coef = TRUE, cor.method = "pearson",
           xlab = "DGCR5 transcript count", ylab = "IFNB transcript count", title = "IFNB vs. DGCR5 COPD infected (18hrs) samples")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text( size=14),
        axis.text.y = element_text( size=14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        plot.title = element_text(size = (14), hjust=0.5))

ggsave(filename = "../../../figs/lncRNAs/DGCR5_IFNB_correlate_COPD_virus_18_samples.pdf",
       plot = last_plot(),
       device = "pdf",
       height = 8,
       width = 10,
       scale = 1,
       dpi = 400
)

DGCR5_IFNB_NHB_virus_18_df <- DGCR5_IFNB_df[NHB_virus_18_idx,]
ggscatter( DGCR5_IFNB_NHB_virus_18_df, x ="DGCR5", y = "IFNB", 
           add = "reg.line", conf.int = TRUE, 
           cor.coef = TRUE, cor.method = "pearson",
           xlab = "DGCR5 transcript count", ylab = "IFNB transcript count", title = "IFNB vs. DGCR5 NHB infected (18hrs) samples")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text( size=14),
        axis.text.y = element_text( size=14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        plot.title = element_text(size = (14), hjust=0.5))
ggsave(filename = "../../../figs/lncRNAs/DGCR5_IFNB_correlate_NHB_virus_18_samples.pdf",
       plot = last_plot(),
       device = "pdf",
       height = 8,
       width = 10,
       scale = 1,
       dpi = 400
)
# DGCR5_IFNB_COPD_virus_18_df <- data.frame(matrix(nrow = length(COPD_virus_18_idx)))
# DGCR5_IFNB_NHB_virus_18_df <- data.frame(matrix(nrow = length(NHB_virus_18_idx)))

DGCR5_IFNB_df$DGCR5 <- as.matrix(DGCR5_counts)[1,]
DGCR5_IFNB_df$IFNB <- as.matrix(IFNB_counts)[1,]

sample_data_copd_nhb$donor_ID[c(33,37,75, 80,114,115,118,119,152,153,154,156,157,159,161,162)]


##--------------------------------------- Identifying virus-specific biomarkers

# viral specific genes
fold_thr <- 1
pval_thr <- 0.05

# 4hrs is not interesting, unless 
virus_sig_nhb_4hrs <- with(res_nhb_v4, which(padj < pval_thr & abs(log2FoldChange) > fold_thr)) # there was a typo here before -- was previously not in absolute value
virus_sig_nhb_18hrs <- with(res_nhb_v18, which(padj < pval_thr & abs(log2FoldChange) > fold_thr)) # there was a typo here before -- was previously not in absolute value

virus_sig_copd_4hrs <- with(res_copd_v4, which(padj < pval_thr & abs(log2FoldChange) > fold_thr)) # there was a typo here before -- was previously not in absolute value
virus_sig_copd_18hrs <- with(res_copd_v18, which(padj < pval_thr & abs(log2FoldChange) > fold_thr)) # there was a typo here before -- was previously not in absolute value

copd_nhb_diff_sig_4hrs <- with(res_copd_nhb_contr4, which(padj < pval_thr & abs(log2FoldChange) > fold_thr)) # there was a typo here before -- was previously not in absolute value
copd_nhb_diff_sig18hrs <- with(res_copd_nhb_contr18, which(padj < pval_thr & abs(log2FoldChange) > fold_thr)) # there was a typo here before -- was previously not in absolute value

poly_ic_sig_copd_18hrs <- with(res_copd_poly_ic18, which(padj < pval_thr & abs(log2FoldChange) > fold_thr)) # there was a typo here before -- was previously not in absolute value
poly_ic_sig_nhb_18hrs <- with(res_nhb_poly_ic18, which(padj < pval_thr & abs(log2FoldChange) > fold_thr)) # there was a typo here before -- was previously not in absolute value
write_csv(data.frame(res_nhb_poly_ic18),"test.csv")
# extracting genes that are differentially expressed in influenza, but not poly-IC challenged models
virus_specific_nhb_18hrs <- setdiff(virus_sig_nhb_18hrs, poly_ic_sig_nhb_18hrs)
virus_specific_copd_18hrs <- setdiff(virus_sig_copd_18hrs, poly_ic_sig_copd_18hrs)

virus_specific_copd_v_nhb_18hrs <- setdiff(virus_specific_copd_18hrs,virus_specific_nhb_18hrs)
virus_specific_nhb_v_copd_18hrs <- setdiff(virus_specific_nhb_18hrs, virus_specific_copd_18hrs)

virus_specific_nhb_genes_18hrs <- gene_ids_copd_nhb[virus_specific_nhb_18hrs]
virus_specific_copd_genes_18hrs <- gene_ids_copd_nhb[virus_specific_copd_18hrs]

virus_nhb_genes <- gene_ids_copd_nhb[virus_sig_nhb_18hrs]
virus_copd_genes <-  gene_ids_copd_nhb[virus_sig_copd_18hrs]
virus_copd_nhb_both <- union(virus_nhb_genes,virus_copd_genes)
polyIC_copd_genes_18hrs <- gene_ids_copd_nhb[poly_ic_sig_copd_18hrs]
polyIC_nhb_genes_18hrs <- gene_ids_copd_nhb[poly_ic_sig_nhb_18hrs]

copd_nhb_genes_4hrs <- gene_ids_copd_nhb[copd_nhb_diff_sig_4hrs]

copd_specific_virus_genes_18hrs <- gene_ids_copd_nhb[virus_specific_copd_v_nhb_18hrs]
nhb_specific_virus_genes_18hrs <- gene_ids_copd_nhb[virus_specific_nhb_v_copd_18hrs]
copd_diff_genes <- gene_ids_copd_nhb[copd_nhb_diff_sig18hrs]
test2 <- intersect(copd_diff_genes,virus_copd_nhb_both)
#################### Volcano plots

source("util/diff_expression_pathway_analysis_functions.R") # script for gene- and pathway-level analysis
copd_specific_and_shared <- unique(union(setdiff(virus_copd_genes, virus_nhb_genes), sample(intersect(virus_copd_genes,virus_nhb_genes),20)))
nhb_specific_and_shared <- unique(union(setdiff(virus_nhb_genes, virus_copd_genes), sample(intersect(virus_copd_genes,virus_nhb_genes),30)))

volcano_nhb_copd_output_dir = "../../../figs/COPD_NHB_RNAseq/differential_expression/"
plot_volcano(res_copd_nhb_contr4, -2, 7, pcutoff = 0.05, paste(volcano_nhb_copd_output_dir, 'volcano_plot_copd_nhb_contr4',sep=''))
plot_volcano(res_copd_nhb_contr18, -2, 7, pcutoff = 0.05, paste(volcano_nhb_copd_output_dir, 'volcano_plot_copd_nhb_contr18',sep=''))

plot_volcano(res_nhb_v4, -2, 7, pcutoff = 0.05, paste(volcano_nhb_copd_output_dir, 'volcano_plot_nhb_v4',sep=''))
plot_volcano (res_nhb_v18, -2, 7, pcutoff = 0.05, paste(volcano_nhb_copd_output_dir, 'volcano_plot_nhb_v18',sep=''))

plot_volcano(res_copd_v4, -2, 7, pcutoff = 0.05, paste(volcano_nhb_copd_output_dir, 'volcano_plot_copd_v4',sep=''))
plot_volcano(res_copd_v18, -2, 7, pcutoff = 0.05, paste(volcano_nhb_copd_output_dir, 'volcano_plot_copd_v18',sep=''))

plot_volcano(res_copd_nhb_v4, -2, 7, pcutoff = 0.05, paste(volcano_nhb_copd_output_dir, 'volcano_plot_nhb_copd_v4',sep=''))
plot_volcano(res_copd_nhb_v18, -2, 7, pcutoff = 0.05, paste(volcano_nhb_copd_output_dir, 'volcano_plot_nhb_copd_v18',sep=''))

plot_volcano(res_copd_poly_ic18, -2, 7, pcutoff = 0.05, paste(volcano_nhb_copd_output_dir, 'volcano_plot_copd_poly_ic_18',sep=''))
plot_volcano(res_nhb_poly_ic18, -2, 7, pcutoff = 0.05, paste(volcano_nhb_copd_output_dir, 'volcano_plot_nhb_poly_ic_18',sep=''))

plot_volcano(res_copd_nhb_v18, -2, 7, pcutoff = 0.05, 'test')

# res_copd_nhb_v18[which(res_copd_nhb_v18$gene == 'DGCR5'), ]
# 
# write.csv(res_copd_nhb_contr18, file = "test.csv")
write.csv(res_copd_nhb_contr4$gene[copd_nhb_diff_sig_4hrs], file = "COPD_NHB_DEG.csv")

# write.table(virus_specific_genes, file="../data/gene_lists/virus_specific_biomarkers_chip.txt", sep="\t", col.names = F, row.names = F, quote=F)

##--------------------------------------- Tissue specificty with Stouffer coefficient

# biogps data ---- tissue and cell expression data
biogps_data = load_bioGPS()

avg_tissue_exp <- compute_avg_tissue_exp(biogps_data)

output_dir_stouffer <- "../../../figs/COPD_NHB_RNAseq/tissue_specificity/"
lung_flu_nhb_biomarkers_list_18hrs <- lung_specific_biomarkers(avg_tissue_exp, virus_nhb_genes,biogps_data, stouf_thr=3, num_tissues=2, output_dir_stouffer)
lung_flu_nhb_biomarkers_18hrs <- lung_flu_nhb_biomarkers_list_18hrs$lung_spcific
lung_nhb_genes_18hrs <- lung_flu_nhb_biomarkers_list_18hrs$lung_genes

lung_flu_copd_biomarkers_list_18hrs <- lung_specific_biomarkers(avg_tissue_exp, virus_copd_genes,biogps_data, stouf_thr=3, num_tissues=2, output_dir_stouffer)
lung_flu_copd_biomarkers_18hrs <- lung_flu_copd_biomarkers_list_18hrs$lung_spcific
lung_copd_genes_18hrs <- lung_flu_copd_biomarkers_list_18hrs$lung_genes

# lung_flu_copd_specific_biomarkers_list_18hrs <- lung_specific_biomarkers(avg_tissue_exp, copd_specific_virus_genes_18hrs,biogps_data, stouf_thr=3, num_tissues=2, output_dir_stouffer)
# lung_flu_copd_specificbiomarkers_18hrs <- lung_flu_copd_specific_biomarkers_list_18hrs$lung_spcific
# lung_copd_specific_genes_18hrs <- lung_flu_copd_specific_biomarkers_list_18hrs$lung_genes
# write.table(lung_flu_biomarkers, file="../data/analysis/lung_specific_flu_biomarkers.txt", sep="\t", col.names = F, row.names = F, quote=F)

# ##--------------------------------------- Gene set enrichment analysis (GAGE)
#count_data_copd_nhb
count_data_copd_nhb_new_colNames <- count_data_copd_nhb
colnames(count_data_copd_nhb_new_colNames) <- sample_data_copd_nhb$condition_group_time
# counts_copd_nhb_normalized
cn <- colnames(count_data_copd_nhb_new_colNames)
copd_contr_4_idx <- which(cn %in% "COPD_control_4")
copd_virus_4_idx <- which(cn %in% "COPD_virus_4")

nhb_contr_4_idx <- which(cn %in% "NHB_control_4")
nhb_virus_4_idx <- which(cn %in% "NHB_virus_4")

copd_contr_18_idx <- which(cn %in% "COPD_control_18")
copd_virus_18_idx <- which(cn %in% "COPD_virus_18")

nhb_contr_18_idx <- which(cn %in% "NHB_control_18")
nhb_virus_18_idx <- which(cn %in% "NHB_virus_18")

copd_poly_ic_18_idx <- which(cn %in% "COPD_poly_ic_18")
nhb_poly_ic_18_idx <- which(cn %in% "NHB_poly_ic_18")

rownames(counts_chip_normalized) <- gkey_all$ENTREZID
# GSEA_Gage(count_data_chip_df, control_idx_v18, case_idx_v18, gene_set_kegg = kegg.gs, output_loc="../../../figs/original_RNA-seq/pathway_analysis/GSEA_GAGE_chip_original", time_point = "18hrs", num_paths = 10, draw_map = FALSE)
# GSEA_Gage(count_data_chip_df, control_idx_v48, case_idx_v48, gene_set_kegg = kegg.gs, output_loc="../../../figs/original_RNA-seq/pathway_analysis/GSEA_GAGE_chip_original", time_point = "48hrs", num_paths = 10, draw_map = FALSE)
# GSEA_Gage(count_data_chip_df, control_idx_IC, case_idx_IC, gene_set_kegg = kegg.gs, output_loc="../../../figs/original_RNA-seq/pathway_analysis/GSEA_GAGE_chip_original_polyIC", time_point = "", num_paths = 10, draw_map = FALSE)

# res_copd_nhb_contr4
# res_copd_nhb_contr18
# 
# res_copd_v4
# res_copd_v18
# 
# res_nhb_v4
# res_nhb_v18
# 
# res_copd_nhb_v4
# res_copd_nhb_v18
# 
# res_copd_poly_ic18
# res_nhb_poly_ic18

source("util/diff_expression_pathway_analysis_functions.R") # script for gene- and pathway-level analysis

rownames(counts_copd_nhb_normalized) <- gkey_all$ENTREZID
# GSEA_Gage(count_data_chip_df, control_idx_v18, case_idx_v18, gene_set_kegg = kegg.gs, output_loc="../../../figs/original_RNA-seq/pathway_analysis/GSEA_GAGE_chip_original", time_point = "18hrs", num_paths = 10, draw_map = FALSE)
res_nhb_v18_entrez <-res_nhb_v18
res_nhb_v18_entrez$entrez = gkey_all$ENTREZID

res_copd_v18_entrez <-res_copd_v18
res_copd_v18_entrez$entrez = gkey_all$ENTREZID

res_copd_v_nhb_4_entrez <-res_copd_nhb_v4
res_copd_v_nhb_4_entrez$entrez = gkey_all$ENTREZID

res_copd_v_nhb_18_entrez <- res_copd_nhb_v18
res_copd_v_nhb_18_entrez$entrez = gkey_all$ENTREZID

res_copd_nhb_contr18_entrez <- res_copd_nhb_contr18
res_copd_nhb_contr18_entrez$entrez <- gkey_all$ENTREZID


enriched_paths_nhb_v18hrs <-GSEA_Gage_old(res_nhb_v18_entrez, "DESeq2", G = kegg.sets.hs, "../../../figs/COPD_NHB_RNAseq/pathway_analysis/GSEA_GAGE_chip_nhb_v", "18hrs", draw_map = FALSE)
enriched_paths_copd_v18hrs <-GSEA_Gage_old(res_copd_v18_entrez, "DESeq2", G = kegg.sets.hs, "../../../figs/COPD_NHB_RNAseq/pathway_analysis/GSEA_GAGE_chip_copd_v", "18hrs", draw_map = FALSE)
enriched_paths_copd_vs_nhb_18hrs <-GSEA_Gage_old(res_copd_v_nhb_18_entrez, "DESeq2", G = kegg.sets.hs, "../../../figs/COPD_NHB_RNAseq/pathway_analysis/GSEA_GAGE_chip_copd_vs_nhb", "18hrs", draw_map = FALSE)
enriched_paths_copd_vs_nhb_v18hrs <-GSEA_Gage_old(res_copd_nhb_v18_entrez, "DESeq2", G = kegg.sets.hs, "../../../figs/COPD_NHB_RNAseq/pathway_analysis/GSEA_GAGE_chip_copd_vs_nhb", "18hrs", draw_map = FALSE)

enriched_paths_res_copd_nhb_contr18_entrez <-GSEA_Gage_old(res_copd_nhb_contr18_entrez, "DESeq2", G = kegg.sets.hs, "test", "18hrs", draw_map = FALSE, num_paths = 10)

# enriched_paths_nhb_v4hrs <- GSEA_Gage(counts_copd_nhb_normalized, nhb_contr_18_idx, nhb_virus_18_idx, gene_set_kegg = kegg.gs, output_loc="../../../figs/COPD_NHB_RNAseq/pathway_analysis/GSEA_GAGE_chip_nhb_virus", time_point = "18hrs", num_paths = 3, draw_map = FALSE)

# enriched_paths_chip_48hrs <- GSEA_Gage(counts_chip_normalized, control_idx_v48, case_idx_v48, gene_set_kegg = kegg.gs, output_loc="../../../figs/original_RNA-seq/pathway_analysis/GSEA_GAGE_chip_original", time_point = "48hrs", num_paths = 12, draw_map = FALSE)
# enriched_paths_chip_polyic <- GSEA_Gage(counts_chip_normalized, control_idx_IC, case_idx_IC, gene_set_kegg = kegg.gs, output_loc="../../../figs/original_RNA-seq/pathway_analysis/GSEA_GAGE_chip_original_polyIC", time_point = "", num_paths = 12, draw_map = FALSE)


##--------------------------------------- Zaas et al differential expression analysis with limma (for microarray expression data)

# DREAM respiratory challenge data ----
zaas_data <-load_zaas()
E <- zaas_data$expression    # expression matrix
G <- zaas_data$genes # gene data
S <- zaas_data$metadata # sample metadata

DE_Zaas <- DE_limma_Zaas(S, E)
xx1 <- DE_Zaas$flu
xx2 <- DE_Zaas$rsv
xx3 <- DE_Zaas$hrv

diff_genes_flu <- which(xx1$adj.P.Val < 0.05 & abs(xx1$logFC) > 1)
diff_genes_rsv <- which(xx2$adj.P.Val < 0.05 & abs(xx2$logFC) > 1)
diff_genes_hrv <- which(xx3$adj.P.Val < 0.05 & abs(xx3$logFC) > 1)
diff_genes_flu_all <- diff_genes_flu

# flu biomarkers ----
# these are genes that are differentially expressed in flu but not
# in the other 2 infection models
diff_genes_flu <- setdiff(diff_genes_flu,
                          union(diff_genes_hrv, diff_genes_rsv))

diff_res_flu_old <- tibble::tibble(gene_id=G$ENTREZID,
                                   gene_symbol=G$SYMBOL,
                                   fold_change=xx1$logFC,
                                   fdr_pvalue=xx1$adj.P.Val)

diff_res_flu <- tibble::tibble(gene_id=G$ENTREZID[diff_genes_flu],
                               gene_symbol=G$SYMBOL[diff_genes_flu],
                               fold_change=xx1$logFC[diff_genes_flu],
                               fdr_pvalue=xx1$adj.P.Val[diff_genes_flu])

diff_res_flu_all <- tibble::tibble(gene_id=G$ENTREZID[diff_genes_flu_all],
                                   gene_symbol=G$SYMBOL[diff_genes_flu_all],
                                   fold_change=xx1$logFC[diff_genes_flu_all],
                                   fdr_pvalue=xx1$adj.P.Val[diff_genes_flu_all])

# zaas all  flu 
# res_zaas_all_gene_annotation <- diff_res_flu_old
# diff_res_flu_all_all_annots_to_add <- gkey_all$GENENAME[match(diff_res_flu_old$gene_symbol, gkey_all$SYMBOL)]
# res_zaas_all_gene_annotation["gene_annotation"] <- diff_res_flu_all_all_annots_to_add
# write.csv(res_zaas_all_gene_annotation, file = "../../../data/RNA-seq_data_chip/original_chip_data/DE_comparisons/res_zaas_all_genes_annotation.csv")
# 
# # zaas all DE flu 
# res_zaas_all_DE_gene_annotation <- diff_res_flu_all
# diff_res_flu_all_annots_to_add <- gkey_all$GENENAME[match(diff_res_flu_all$gene_symbol, gkey_all$SYMBOL)]
# res_zaas_all_DE_gene_annotation["gene_annotation"] <- diff_res_flu_all_annots_to_add
# # write.csv(res_zaas_all_DE_gene_annotation, file = "../../../data/RNA-seq_data_chip/original_chip_data/DE_comparisons/res_zaas_all_annotations.csv")
# 
# # zaas specific 
# res_zaas_specific_DE_gene_annotation <- diff_res_flu
# diff_res_flu_annots_to_add <- gkey_all$GENENAME[match(diff_res_flu$gene_symbol, gkey_all$SYMBOL)]
# res_zaas_specific_DE_gene_annotation["gene_annotation"] <- diff_res_flu_annots_to_add
# # write.csv(res_zaas_specific_DE_gene_annotation, file = "../../../data/RNA-seq_data_chip/original_chip_data/DE_comparisons/res_zaas_specific_annotations.csv")
#########
enriched_paths_flu <- GSEA_Gage_old(xx1,DEmethod = "limma", G, time_point="","../../../figs/COPD_NHB_RNAseq/pathway_analysis/GSEA_GAGE_Zaas_flu")
# enriched_paths_rsv <- GSEA_Gage_old(xx2,DEmethod = "limma", G, "../../../figs/original_RNA-seq/pathway_analysis/GSEA_GAGE_Zaas_rsv")
# enriched_paths_hsv <- GSEA_Gage_old(xx3,DEmethod = "limma", G, "../../../figs/original_RNA-seq/pathway_analysis/GSEA_GAGE_Zaas_hsv")

######
#--------------------------------------- comparing pathways (GAGE)
output_dir_compare_zaas_chip_copd = "../../../figs/COPD_NHB_RNAseq/pathway_analysis/"
# compare_chip_Zaas(enriched_paths_chip_18hrs, enriched_paths_flu, direction="upregulated", num_paths=15,out_name_compare_18hrs)
# compare_chip_Zaas(enriched_paths_chip_48hrs, enriched_paths_flu, direction="upregulated", num_paths=15,out_name_compare_48hrs)

source("util/diff_expression_pathway_analysis_functions.R") # script for gene- and pathway-level analysis

compare_chip_copd_Zaas(enriched_paths_nhb_v18hrs, enriched_paths_copd_v18hrs, enriched_paths_flu, direction="downregulated", num_paths=10,output_dir_compare_zaas_chip_copd)

####
## ---------------------------- Lung-specific genes
# plot ----
source("util/diff_expression_pathway_analysis_functions.R") # script for gene- and pathway-level analysis

b2_18 <- res_nhb_v18$log2FoldChange[which(gene_ids_copd_nhb %in% lung_flu_nhb_biomarkers_18hrs)] # <-- chip nhb data
b3_18 <- res_copd_v18$log2FoldChange[which(gene_ids_copd_nhb %in% lung_flu_copd_biomarkers_18hrs)] # <-- chip nhb data

#b3 <- y1$logFC[which(gene_mappings_copd$SYMBOL %in% b1)] # <-- copd data
b4_18 <- xx1$logFC[which(G$SYMBOL %in% lung_flu_copd_biomarkers_18hrs)] # <-- zaas data, flu
b5_18 <- xx2$logFC[which(G$SYMBOL %in% lung_flu_copd_biomarkers_18hrs)] # <-- zaas data, rsv
b6_18 <- xx3$logFC[which(G$SYMBOL %in% lung_flu_copd_biomarkers_18hrs)] # <-- zaas data, hrv

dge_groups <- tibble(group = c(rep("this study (NHB)", length(b2_18)),
                               rep("this study (COPD)", length(b3_18)),
                               #rep("copd data (benam et al)", length(b3)),
                               rep("flu (zaas et al. 2009)", length(b4_18)),
                               rep("rsv (zaas et al. 2009)", length(b5_18)),
                               rep("hrv (zaas et al. 2009)", length(b6_18))),
                     gene_name = c(gene_ids_copd_nhb[gene_ids_copd_nhb %in% lung_flu_copd_biomarkers_18hrs], 
                                   gene_ids_copd_nhb[gene_ids_copd_nhb %in% lung_flu_copd_biomarkers_18hrs],
                                   #gene_mappings_copd$SYMBOL[which(gene_mappings_copd$SYMBOL %in% b1)],
                                   G$SYMBOL[G$SYMBOL %in% lung_flu_copd_biomarkers_18hrs],
                                   G$SYMBOL[G$SYMBOL %in% lung_flu_copd_biomarkers_18hrs],
                                   G$SYMBOL[G$SYMBOL %in% lung_flu_copd_biomarkers_18hrs]),
                     fold_data = c(b2_18, b3_18,b4_18, b5_18, b6_18))

dge_groups %>% 
  dplyr::filter(., group != "copd data (benam et al)", fold_data > 0.90) %>%
  ggplot() + 
  geom_point(aes(x = group, y = fold_data, color = group), size = 7, alpha = 1) + 
  scale_color_manual(values = c('deeppink', 'darkred','green', 'blue', 'orange')) +
  facet_wrap(. ~ gene_name) + 
  labs(x = NULL, y = "log2(fold difference to control)") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 12, color = "black", angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        legend.position = "none")

ggsave(filename = "../../../figs/COPD_NHB_RNAseq/tissue_specificity/lung_specific_biomarkers_18hrs.pdf",
       plot = last_plot(), 
       device = "pdf", 
       height = 8,
       width = 10,
       scale = 1, 
       dpi = 400)

ggsave(filename = "../../../figs/COPD_NHB_RNAseq/tissue_specificity/lung_specific_biomarkers_18hrs.png",
       plot = last_plot(), 
       device = "png", 
       height = 8,
       width = 10,
       scale = 1, 
       dpi = 400)
##--------------------------------------- Signal peptides / secreted biomarkers

# signal peptides ----
sig_pep <- read_xlsx("../../../data/secretome/2019-03-07_signal_peptides_human.xlsx")

sig_peps_names <- AnnotationDbi::select(x = org.Hs.eg.db, keys = as.character(sig_pep$`Accession Number`), 
                                        keytype = "UNIPROT", columns = c("ENTREZID", "SYMBOL"))
sig_peps_names <- sig_peps_names[which(!is.na(sig_peps_names$ENTREZID)), ]

genes_to_compare <- union(virus_nhb_genes,virus_nhb_genes)
# differentially abundant signal peptides in our data
virus_sig_peps_nhb_18 <- intersect(sig_peps_names$SYMBOL, genes_to_compare)
virus_sig_peps_copd_18 <- intersect(sig_peps_names$SYMBOL, genes_to_compare)

# secreted proteins ----
# based on HPA predicted data
hpa_pred <- read_delim("../../../data/secretome/2019-03-07_HPA_predicted_secreted_proteins.tsv", delim = "\t")
hpa_sec <- hpa_pred[grep("secreted", hpa_pred$`Protein class`), ]

# which virus secreted proteins are differentially abundant 
virus_secp_nhb_18 <- intersect(hpa_pred$Gene, genes_to_compare)
virus_secp_copd_18 <- intersect(hpa_pred$Gene, genes_to_compare)
# predicted secreted biomarkers ----
virus_biom_sec_nhb_18 <- union(virus_sig_peps_nhb_18, virus_secp_nhb_18)
virus_biom_sec_copd_18 <- union(virus_sig_peps_copd_18, virus_secp_copd_18)

length(virus_biom_sec_nhb_18)
length(virus_biom_sec_copd_18)

# virus_biom_sec # proteins that are predicted to be secreted and have signal peptides

# write.table(virus_biom_sec, file="../data/gene_lists/virus_secreted_signal_chip.txt", sep="\t", col.names = F, row.names = F, quote=F)
# write.table(virus_secp_48hrs, file="../../data/gene_lists/virus_secreted_signal_chip_48hrs.txt", sep="\t", col.names = F, row.names = F, quote=F)

virus_biom_sec_18 <- union(virus_biom_sec_nhb_18, virus_biom_sec_copd_18)

# secreted virus specific ----
d2_18 <- res_nhb_v18$log2FoldChange[which(gene_ids_copd_nhb %in% virus_biom_sec_18)] # <-- chip data
d3_18 <- res_copd_v18$log2FoldChange[which(gene_ids_copd_nhb %in% virus_biom_sec_18)] # <-- chip data

d4_18<- xx1$logFC[which(G$SYMBOL %in% virus_biom_sec_18)] # <-- zaas data, flu
d5_18 <- xx2$logFC[which(G$SYMBOL %in% virus_biom_sec_18)] # <-- zaas data, rsv
d6_18 <- xx3$logFC[which(G$SYMBOL %in% virus_biom_sec_18)] # <-- zaas data, hrv

sec_genes_18 <- tibble(group = c(rep("this study (NHB patients)", length(d2_18)),
                                 rep("this study (COPD patients)", length(d2_18)),
                                 rep("flu infection data (zaas et al)", length(d4_18)),
                                 rep("rsv infection data (zaas et al)", length(d5_18)),
                                 rep("hrv infection data (zaas et al)", length(d6_18))),
                       gene_name = c(gene_ids_copd_nhb[which(gene_ids_copd_nhb %in% virus_biom_sec_18)], 
                                     gene_ids_copd_nhb[which(gene_ids_copd_nhb %in% virus_biom_sec_18)], 
                                     G$SYMBOL[which(G$SYMBOL %in% virus_biom_sec_18)],
                                     G$SYMBOL[which(G$SYMBOL %in% virus_biom_sec_18)],
                                     G$SYMBOL[which(G$SYMBOL %in% virus_biom_sec_18)]),
                       fold_data = c(d2_18,d3_18, d4_18, d5_18, d6_18))

sec_genes_18 %>% 
  ggplot() + 
  geom_point(aes(x = forcats::fct_reorder(gene_name, fold_data, .desc = TRUE), y = fold_data, color = group), size =7, alpha = 1.0) + 
  geom_hline(yintercept = 1, color = "black", lty = 2) +geom_hline(yintercept = -1, color = "black", lty = 2)+
  # scale_color_manual(values = c("purple", "orange")) +
  labs(x = "Secreted biomarker", y = "log2(fold difference to control)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "bottom")+ 
  scale_color_manual(values = c('deeppink', 'darkred','green', 'blue', 'orange'))

ggsave(filename = "../../../figs/COPD_NHB_RNAseq/differential_expression/secreted_biomarkers_18hrs.pdf",
       plot = last_plot(), 
       device = "pdf", 
       height = 8,
       width = 10,
       scale = 1, 
       dpi = 400)

ggsave(filename = "../../../figs/COPD_NHB_RNAseq/differential_expression/secreted_biomarkers_18hrs.png",
       plot = last_plot(), 
       device = "png", 
       height = 8,
       width = 10,
       scale = 1, 
       dpi = 400)

# secreted virus specific ----
d2_48 <- res_v48$log2FoldChange[which(gene_ids_chip %in% virus_biom_sec_48)] # <-- chip data
d4_48 <- xx1$logFC[which(G$SYMBOL %in% virus_biom_sec_48)] # <-- zaas data, flu
d5_48 <- xx2$logFC[which(G$SYMBOL %in% virus_biom_sec_48)] # <-- zaas data, rsv
d6_48 <- xx3$logFC[which(G$SYMBOL %in% virus_biom_sec_48)] # <-- zaas data, hrv

sec_genes_48 <- tibble(group = c(rep("this study", length(d2_48)),
                                 rep("flu infection data (zaas et al)", length(d4_48)),
                                 rep("rsv infection data (zaas et al)", length(d5_48)),
                                 rep("hrv infection data (zaas et al)", length(d6_48))),
                       gene_name = c(gene_ids_chip[which(gene_ids_chip %in% virus_biom_sec_48)], 
                                     G$SYMBOL[which(G$SYMBOL %in% virus_biom_sec_48)],
                                     G$SYMBOL[which(G$SYMBOL %in% virus_biom_sec_48)],
                                     G$SYMBOL[which(G$SYMBOL %in% virus_biom_sec_48)]),
                       fold_data = c(d2_48, d4_48, d5_48, d6_48))

sec_genes_48 %>% 
  ggplot() + 
  geom_point(aes(x = forcats::fct_reorder(gene_name, fold_data, .desc = TRUE), y = fold_data, color = group), size = 3, alpha = 0.5) + 
  geom_hline(yintercept = 1, color = "black", lty = 2) + geom_hline(yintercept = -1, color = "black", lty = 2)+ 
  scale_color_manual(values = c('deeppink', 'darkred','green', 'blue', 'orange'))
  # scale_color_manual(values = c("purple", "orange")) +
  labs(x = "Secreted biomarker", y = "log2(fold difference to control)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 20),
        axis.text.y = element_text(size = 12),
        legend.position = "bottom")

ggsave(filename = "../../../figs/original_RNA-seq/differential_expression/secreted_biomarkers_48hrs.pdf",
       plot = last_plot(), 
       device = "pdf", 
       height = 8,
       width = 10,
       scale = 1, 
       dpi = 400)

ggsave(filename = "../../../figs/original_RNA-seq/differential_expression/secreted_biomarkers_48hrs.png",
       plot = last_plot(), 
       device = "png", 
       height = 8,
       width = 10,
       scale = 1, 
       dpi = 400)

res_v18_sec <- res_v18[which(gene_ids_chip %in% virus_biom_sec_18),] %>% 
  tibble::add_column(., annotation = subset(gkey_all[,3],  gkey_all[,1] %in% virus_biom_sec_18))
res_v48_sec <- res_v48[which(gene_ids_chip %in% virus_biom_sec_48),] %>% 
  tibble::add_column(., annotation = subset(gkey_all[,3],  gkey_all[,1] %in% virus_biom_sec_48))
# ggsave(filename = "../figs/secreted_biomarkers.svg",
#        plot = last_plot(), 
#        device = "svg", 
#        height = 8,
#        width = 10,
#        scale = 1, 
#        dpi = 400)
# 
# write.csv(res_v18_sec, file = "../../../data/RNA-seq_data_chip/original_chip_data/secreted_biomarkers/res_v18_sec.csv")
# write.csv(res_v48_sec, file = "../../../data/RNA-seq_data_chip/original_chip_data/secreted_biomarkers/res_v48_sec.csv")
