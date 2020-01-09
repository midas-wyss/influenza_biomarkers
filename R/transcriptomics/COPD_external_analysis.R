library("org.Hs.eg.db") # remember to install it if you don't have it already
library("dplyr")
library(tidyverse) # for organizing data structures
library(DESeq2) # differential expression (RNA-seq)
library(EnhancedVolcano)

### ------------------ from (https://respiratory-research.biomedcentral.com/articles/10.1186/s12931-019-1032-z)
# Essentially RNA-seq on three different tissue types in COPD patients
source("util/load_transcriptomics_data_functions.R")
source("util/preprocess_transcriptomics_data_functions.R")
source("util/diff_expression_pathway_analysis_functions.R")

COPD_external_data <- load_COPD_external()
copd_count_data_original <- COPD_external_data$copd_count_data
copd_genes_ENSEMBL <- COPD_external_data$genes
copd_metadata <- COPD_external_data$metadata
# COPD_gene_annots <- AnnotationDbi::select(org.Hs.eg.db, keys = COPD_gene_ids_EMSEMBL, keytype = "ENSEMBL", columns=c("SYMBOL")) # , "ENTREZID", "GENENAME"
copd_external_data_processed <- preprocess_copd_external_data(copd_count_data_original, copd_genes_ENSEMBL, org.Hs.eg.db)

copd_count_processed <- copd_external_data_processed$count_data
copd_genes_processed <- copd_external_data_processed$gene_ids
copd_filtered_genes <- copd_external_data_processed$genes_filtered
copd_filtered_gkey <- copd_external_data_processed$gkey_filtered

copd_gene_symbols <- copd_filtered_gkey$SYMBOL
dds_chip <- DESeqDataSetFromMatrix(countData = copd_count_processed,
                                   colData = copd_metadata,
                                   design = ~ group_location)
dds_chip <- DESeq(object = dds_chip)


contrast_alveolar_macrophage <- c("group_location", "copd_MO", "control_MO")
contrast_bronchial_epithelial <- c("group_location", "copd_BE", "control_BE")
contrast_whole_blood <- c("group_location", "copd_WB", "control_WB")

res_alveolar_macrophage <- compare_group_DESeq2(dds_chip, contrast_alveolar_macrophage, copd_gene_symbols, MHT = "fdr")
res_bronchial_epithelial <- compare_group_DESeq2(dds_chip, contrast_bronchial_epithelial, copd_gene_symbols, MHT = "fdr")
res_whole_blood  <- compare_group_DESeq2(dds_chip, contrast_whole_blood, copd_gene_symbols, MHT = "fdr")


DE_MO <- with(res_alveolar_macrophage, which(padj < 0.05 & abs(log2FoldChange) > 1)) # there was a typo here before -- was previously not in absolute value
DE_BE <- with(res_bronchial_epithelial, which(padj < 0.05 & abs(log2FoldChange) > 1)) # there was a typo here before -- was previously not in absolute value
DE_WB <- with(res_whole_blood, which(padj < 0.05 & abs(log2FoldChange) > 1)) # there was a typo here before -- was previously not in absolute value

# res_MO_all_gene_annotation <- res_alveolar_macrophage
# res_MO_all_annots_to_add <- copd_filtered_gkey$GENENAME[match(res_alveolar_macrophage$gene, copd_filtered_gkey$SYMBOL)]
# res_MO_all_gene_annotation["gene_annotation"] <- diff_MO_all_annots_to_add
# write.csv(res_MO_all_gene_annotation, file = "../../../data/COPD_external_RNAseq/DE_comparisons/res_MO_all_gene_annotation.csv")
# res_MO_DE_only <- res_alveolar_macrophage[DE_MO,]
# res_MO_all_annots_to_add_DE_only <- copd_filtered_gkey$GENENAME[match(res_MO_DE_only$gene, copd_filtered_gkey$SYMBOL)]
# res_MO_DE_only["gene_annotation"] <- res_MO_all_annots_to_add_DE_only
# write.csv(res_MO_DE_only, file = "../../../data/COPD_external_RNAseq/DE_comparisons/res_MO_DE_only.csv")

# res_BE_all_gene_annotation <- res_bronchial_epithelial
# res_BE_all_annots_to_add <- copd_filtered_gkey$GENENAME[match(res_bronchial_epithelial$gene, copd_filtered_gkey$SYMBOL)]
# res_BE_all_gene_annotation["gene_annotation"] <- res_BE_all_annots_to_add
# write.csv(res_BE_all_gene_annotation, file = "../../../data/COPD_external_RNAseq/DE_comparisons/res_BE_all_gene_annotation.csv")

# res_WB_all_gene_annotation <- res_whole_blood
# res_WB_all_annots_to_add <- copd_filtered_gkey$GENENAME[match(res_whole_blood$gene, copd_filtered_gkey$SYMBOL)]
# res_WB_all_gene_annotation["gene_annotation"] <- res_WB_all_annots_to_add
# write.csv(res_WB_all_gene_annotation, file = "../../../data/COPD_external_RNAseq/DE_comparisons/res_WB_all_gene_annotation.csv")

plot_volcano(res_alveolar_macrophage, -2, 7,pcutoff = 0.05, "../../../figs/COPD_external/volcano_plot_COPDexternal_MO")
plot_volcano(res_bronchial_epithelial, -2, 7,pcutoff = 0.05, "../../../figs/COPD_external/volcano_plot_COPDexternal_BE")
plot_volcano(res_whole_blood, -2, 7,pcutoff = 0.05, "../../../figs/COPD_external/volcano_plot_COPDexternal_WB")


res_bronchial_epithelial[DE_BE,]

### ------------------ from (https://www.ncbi.nlm.nih.gov/pubmed/28287180)
library(Biobase)
library(GEOquery)

# load series and platform data from GEO

gset <- getGEO("GSE76925", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL10558", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# set parameters and draw the plot

# dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE76925", '/', annotation(gset), " selected samples", sep ='')
boxplot(exprs(gset), boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)
astibble(exprs(gset))

gset@featureData





