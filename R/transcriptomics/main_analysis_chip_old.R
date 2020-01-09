##--------------------------------------- Influenza biomarkers project overview
# The goal of this project is to computationally identify biomkarers that validate the efficacy of an in-vitro "organ-on-chip" model of influenza infection. 
# Obtaining these biomakers will aid in i) evaluating the model's ability to recapitulate disease physiology ii) identifying novel clinical biomarkers that 
# distinguish influenza from other respiratory diseases and iii) aid in the development of therapeutics tageting influenza. These biomarkers will be recovered 
# from RNA-sequencing data capturing gene expression dynamics in a healthy and infected (COPD) model.  This project is being conducted in collaboration with the 
# Ingber lab and will utilize their acquired RNA-seq data alongside publicly available datasets. 

# RNA-seq is of lung epithelial cells

# Author: Miguel A. Alcantar (initial code provided by Diogo M. Camacho)
# Last updated: 09/26/2019

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

chip_data_original <- load_original_chip_data() # loading all data from organ-on-chip experiments

count_data_chip_original <- chip_data_original$count_data # raw counts
gene_ids_chip_original <- chip_data_original$ge # genes
sample_data_chip <- chip_data_original$metadata # metadata
# writing all gene IDs to a text file
# write.table(gene_ids_chip_original, file="../../../data/gene_lists/gene_ids_chip_all.txt", sep="\t", col.names = F, row.names = F, quote=F)

##--------------------------------------- Data preprocessing 
# Preprocessing will involve:
# removing genes with all zeros for all sample
# remove genes with multiple entrez mappings
# removing genes with <5 counts in 75% of samples

# final gene count and ID data matrices
chip_data_processed <- preprocess_chip_data(count_data_chip_original, gene_ids_chip_original, org.Hs.eg.db)
count_data_chip <- chip_data_processed$count_data # preprocessed count data 
gene_ids_chip <- chip_data_processed$gene_ids # genes remaining after preprocessing
genes_filtered <- chip_data_processed$genes_filtered
gkey_all <- AnnotationDbi::select(x = org.Hs.eg.db, 
                                  keys = gene_ids_chip, 
                                  columns = c("ENTREZID","GENENAME"), 
                                  keytype = "SYMBOL")
## code for saving these data to a csv 
# write.csv(count_data_chip, file = "../../../data/RNA-seq_data_chip/original_chip_data/count_data_organ_on_chip_preprocessed.csv")
# write.table(gene_ids_chip, file= "../../../data/gene_lists/gene_ids_chip_preprocessed.txt", sep="\t", col.names = F, row.names = F, quote=F)
# write.csv(sample_data_chip, file = "../../../data/RNA-seq_data_chip/original_chip_data/organ_on_chip_metadata.csv")
##--------------------------------------- Differential expression analysis with DESeq2

# differential expression ----
# use deseq2
dds_chip <- DESeqDataSetFromMatrix(countData = count_data_chip,
                              colData = sample_data_chip,
                              design = ~ group_time)

dds_chip <- DESeq(object = dds_chip)
counts_chip_normalized <- counts(dds_chip, normalized=TRUE) # normalize counts using sizeFactor; useful for GSEA_GAGE

# conditions we want to compare
contrast_v18 <- c("group_time", "virus_18", "control_18")
contrast_v48 <- c("group_time", "virus_48", "control_48")
contrast_ic18 <- c("group_time", "poly_ic_18", "control_18")
contrast_c48_c18 <- c("group_time", "control_48", "control_18") 

res_v18 <- compare_group_DESeq2(dds_chip, contrast_v18, gene_ids_chip, MHT = "fdr")
res_v48 <- compare_group_DESeq2(dds_chip, contrast_v48, gene_ids_chip, MHT = "fdr")
res_ic18 <- compare_group_DESeq2(dds_chip, contrast_ic18, gene_ids_chip, MHT = "fdr")
res_c48_c18 <- compare_group_DESeq2(dds_chip, contrast_c48_c18, gene_ids_chip, MHT = "fdr")

volcano_output_dir = "../../../figs/original_RNA-seq/differential_expression/"
plot_volcano(res_c48_c18, -2, 7, pcutoff = 0.05, paste(volcano_output_dir, 'volcano_plot_c48_c18',sep=''))
plot_volcano(res_ic18, -2, 7, pcutoff = 0.05,paste(volcano_output_dir,'volcano_plot_ic18',sep=''))
plot_volcano(res_v18, -2, 7, pcutoff = 0.05,paste(volcano_output_dir,'volcano_plot_v18',sep=''))
plot_volcano(res_v48, -2, 7,pcutoff = 0.05,paste(volcano_output_dir,'volcano_plot_v48',sep='') )

# write.csv(res_v18, file = "../../../data/RNA-seq_data_chip/original_chip_data/DE_comparisons/res_v18.csv")
# write.csv(res_v48, file = "../../../data/RNA-seq_data_chip/original_chip_data/DE_comparisons/res_v48.csv")
# write.csv(res_ic18, file = "../../../data/RNA-seq_data_chip/original_chip_data/DE_comparisons/res_ic18.csv")
# write.csv(res_c48_c18, file = "../../../data/RNA-seq_data_chip/original_chip_data/DE_comparisons/res_c48_c18.csv")

##--------------------------------------- Identifying virus-specific biomarkers

# viral specific genes
fold_thr <- 1
pval_thr <- 0.05

virus_sig_18hrs <- with(res_v18, which(padj < pval_thr & abs(log2FoldChange) > fold_thr)) # there was a typo here before -- was previously not in absolute value
ic_sig <- with(res_ic18, which(padj < pval_thr & abs(log2FoldChange) > fold_thr))
controls_sig <- with(res_c48_c18, which(padj < pval_thr & abs(log2FoldChange) > fold_thr)) # it appears that a significant number of genes are differentially expressed in the 48 and 18 hour controls
virus_sig_48hrs <- with(res_v48, which(padj < pval_thr & abs(log2FoldChange) > fold_thr)) 

# extracting genes that are differentially expressed in influenza, but not poly-IC challenged models
virus_specific_18hrs <- setdiff(virus_sig_18hrs, ic_sig)
virus_specific_48hrs <- setdiff(virus_sig_48hrs, ic_sig)

virus_specific_genes_18hrs <- gene_ids_chip[virus_specific_18hrs]
virus_specific_genes_48hrs <- gene_ids_chip[virus_specific_48hrs]
polyIC_genes_18hrs <- gene_ids_chip[ic_sig]
# write.table(virus_specific_genes, file="../data/gene_lists/virus_specific_biomarkers_chip.txt", sep="\t", col.names = F, row.names = F, quote=F)

##--------------------------------------- DE with annotations
# virus_18hrs
res_v18_DE_gene_annotation <- res_v18[virus_sig_18hrs,]
res_v18_annots_to_add <- gkey_all$GENENAME[match(res_v18_DE_gene_annotation$gene, gkey_all$SYMBOL)]
res_v18_DE_gene_annotation["gene_annotation"] <- res_v18_annots_to_add
# write.csv(res_v18_DE_gene_annotation, file = "../../../data/RNA-seq_data_chip/original_chip_data/DE_comparisons/res_v18_DE_annotations.csv")

# virus_18hrs
res_v48_DE_gene_annotation <- res_v48[virus_sig_48hrs,]
res_v48_annots_to_add <- gkey_all$GENENAME[match(res_v48_DE_gene_annotation$gene, gkey_all$SYMBOL)]
res_v48_DE_gene_annotation["gene_annotation"] <- res_v48_annots_to_add
# write.csv(res_v48_DE_gene_annotation, file = "../../../data/RNA-seq_data_chip/original_chip_data/DE_comparisons/res_v48_DE_annotations.csv")

# ic18
res_ic18_DE_gene_annotation <- res_ic18[ic_sig,]
res_ic18_annots_to_add <- gkey_all$GENENAME[match(res_ic18_DE_gene_annotation$gene, gkey_all$SYMBOL)]
res_ic18_DE_gene_annotation["gene_annotation"] <- res_ic18_annots_to_add
# write.csv(res_ic18_DE_gene_annotation, file = "../../../data/RNA-seq_data_chip/original_chip_data/DE_comparisons/res_ic18_DE_annotations.csv")

# c48_18
res_contr_DE_gene_annotation <- res_c48_c18[controls_sig,]
res_contr_annots_to_add <- gkey_all$GENENAME[match(res_contr_DE_gene_annotation$gene, gkey_all$SYMBOL)]
res_contr_DE_gene_annotation["gene_annotation"] <- res_contr_annots_to_add
# write.csv(res_contr_DE_gene_annotation, file = "../../../data/RNA-seq_data_chip/original_chip_data/DE_comparisons/res_contr_DE_annotations.csv")

# virus_specific_genes 18hrs
res_v18_specific_DE_gene_annotation <- res_v18[virus_specific_18hrs,]
res_v18_specific_annots_to_add <- gkey_all$GENENAME[match(res_v18_specific_DE_gene_annotation$gene, gkey_all$SYMBOL)]
res_v18_specific_DE_gene_annotation["gene_annotation"] <- res_v18_specific_annots_to_add
# write.csv(res_v18_specific_DE_gene_annotation, file = "../../../data/RNA-seq_data_chip/original_chip_data/DE_comparisons/res_v18_specific_DE_annotations.csv")

# virus_specific_genes 48hrs
res_v48_specific_DE_gene_annotation <- res_v48[virus_specific_48hrs,]
res_v48_specific_annots_to_add <- gkey_all$GENENAME[match(res_v48_specific_DE_gene_annotation$gene, gkey_all$SYMBOL)]
res_v48_specific_DE_gene_annotation["gene_annotation"] <- res_v48_specific_annots_to_add
# write.csv(res_v48_specific_DE_gene_annotation, file = "../../../data/RNA-seq_data_chip/original_chip_data/DE_comparisons/res_v48_specific_DE_annotations.csv")

##--------------------------------------- Gene set enrichment analysis (GAGE)

# adding entrez IDs to our log2-fold change matrix
res_v18_EID <- res_v18 %>% 
  tibble::add_column(., entrez = subset(gkey_all,  gkey_all[,1] %in% res_v18$gene)[,2])

res_v48_EID <- res_v48 %>% 
  tibble::add_column(., entrez = subset(gkey_all,  gkey_all[,1] %in% res_v48$gene)[,2])
# enriched_paths_chip_18hrs <- GSEA_Gage(res_v18_EID, DEmethod = "DESeq2",  output_loc="../../../figs/original_RNA-seq/pathway_analysis/GSEA_GAGE_chip_original", time_point="18hrs")
# enriched_paths_chip_48hrs <- GSEA_Gage(res_v48_EID, DEmethod = "DESeq2",  output_loc="../../../figs/original_RNA-seq/pathway_analysis/GSEA_GAGE_chip_original", time_point="48hrs")

cn <- colnames(count_data_chip)
control_idx_v18 <- which(colnames(count_data_chip) %in% cn[grep("May1$",cn)][grep("^C", cn[grep("May1$",cn)])])
case_idx_v18 <- which(colnames(count_data_chip) %in% cn[grep("May1$",cn)][grep("^V", cn[grep("May1$",cn)])])

control_idx_v48 <- which(colnames(count_data_chip) %in% cn[grep("May2$",cn)][grep("^C", cn[grep("May2$",cn)])])
case_idx_v48 <- which(colnames(count_data_chip) %in% cn[grep("May2$",cn)][grep("^V", cn[grep("May2$",cn)])])

control_idx_IC <- which(colnames(count_data_chip) %in% cn[grep("May2$",cn)][grep("^C", cn[grep("May2$",cn)])])
case_idx_IC <- which(colnames(count_data_chip) %in% cn[grep("IC", cn)])

rownames(counts_chip_normalized) <- gkey_all$ENTREZID
# GSEA_Gage(count_data_chip_df, control_idx_v18, case_idx_v18, gene_set_kegg = kegg.gs, output_loc="../../../figs/original_RNA-seq/pathway_analysis/GSEA_GAGE_chip_original", time_point = "18hrs", num_paths = 10, draw_map = FALSE)
# GSEA_Gage(count_data_chip_df, control_idx_v48, case_idx_v48, gene_set_kegg = kegg.gs, output_loc="../../../figs/original_RNA-seq/pathway_analysis/GSEA_GAGE_chip_original", time_point = "48hrs", num_paths = 10, draw_map = FALSE)
# GSEA_Gage(count_data_chip_df, control_idx_IC, case_idx_IC, gene_set_kegg = kegg.gs, output_loc="../../../figs/original_RNA-seq/pathway_analysis/GSEA_GAGE_chip_original_polyIC", time_point = "", num_paths = 10, draw_map = FALSE)
enriched_paths_chip_18hrs <- GSEA_Gage(counts_chip_normalized, control_idx_v18, case_idx_v18, gene_set_kegg = kegg.gs, output_loc="../../../figs/original_RNA-seq/pathway_analysis/GSEA_GAGE_chip_original", time_point = "18hrs", num_paths = 12, draw_map = FALSE)
enriched_paths_chip_48hrs <- GSEA_Gage(counts_chip_normalized, control_idx_v48, case_idx_v48, gene_set_kegg = kegg.gs, output_loc="../../../figs/original_RNA-seq/pathway_analysis/GSEA_GAGE_chip_original", time_point = "48hrs", num_paths = 12, draw_map = FALSE)
enriched_paths_chip_polyic <- GSEA_Gage(counts_chip_normalized, control_idx_IC, case_idx_IC, gene_set_kegg = kegg.gs, output_loc="../../../figs/original_RNA-seq/pathway_analysis/GSEA_GAGE_chip_original_polyIC", time_point = "", num_paths = 12, draw_map = FALSE)

##--------------------------------------- gene set enrichment analysis (rpegeos)

chip_enr_18hrs <- rpegeos_enrichment(virus_specific_genes_18hrs, res_v18_EID, num_paths=12, org.Hs.eg.db = org.Hs.eg.db, output_loc="../../../figs/original_RNA-seq/pathway_analysis/GSEA_repegeos_chip_18hrs_original")
chip_enr_48hrs <- rpegeos_enrichment(virus_specific_genes_48hrs, res_v48_EID, num_paths=12, org.Hs.eg.db = org.Hs.eg.db, output_loc="../../../figs/original_RNA-seq/pathway_analysis/GSEA_repegeos_chip_48hrs_original")
chip_enr_polyIC <- rpegeos_enrichment(polyIC_genes_18hrs, res_ic18, num_paths=12, org.Hs.eg.db = org.Hs.eg.db, output_loc="../../../figs/original_RNA-seq/pathway_analysis/GSEA_repegeos_chip_polyIC_original")

##--------------------------------------- Tissue specificty with Stouffer coefficient

# biogps data ---- tissue and cell expression data
biogps_data = load_bioGPS()

avg_tissue_exp <- compute_avg_tissue_exp(biogps_data)

output_dir_stouffer <- "../../../figs/original_RNA-seq/tissue_specificity/"
lung_flu_biomarkers_list_18hrs <- lung_specific_biomarkers(avg_tissue_exp, virus_specific_genes_18hrs,biogps_data, stouf_thr=3, num_tissues=2, output_dir_stouffer)
lung_flu_biomarkers_18hrs <- lung_flu_biomarkers_list_18hrs$lung_spcific
lung_genes_18hrs <- lung_flu_biomarkers_list_18hrs$lung_genes

lung_flu_biomarkers_list_48hrs <- lung_specific_biomarkers(avg_tissue_exp, virus_specific_genes_48hrs,biogps_data, stouf_thr=3, num_tissues=2, output_dir_stouffer)
lung_flu_biomarkers_48hrs <- lung_flu_biomarkers_list_48hrs$lung_spcific
lung_genes_48hrs <- lung_flu_biomarkers_list_48hrs$lung_genes
# write.table(lung_flu_biomarkers, file="../data/analysis/lung_specific_flu_biomarkers.txt", sep="\t", col.names = F, row.names = F, quote=F)

# 
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
res_zaas_all_gene_annotation <- diff_res_flu_old
diff_res_flu_all_all_annots_to_add <- gkey_all$GENENAME[match(diff_res_flu_old$gene_symbol, gkey_all$SYMBOL)]
res_zaas_all_gene_annotation["gene_annotation"] <- diff_res_flu_all_all_annots_to_add
write.csv(res_zaas_all_gene_annotation, file = "../../../data/RNA-seq_data_chip/original_chip_data/DE_comparisons/res_zaas_all_genes_annotation.csv")

# zaas all DE flu 
res_zaas_all_DE_gene_annotation <- diff_res_flu_all
diff_res_flu_all_annots_to_add <- gkey_all$GENENAME[match(diff_res_flu_all$gene_symbol, gkey_all$SYMBOL)]
res_zaas_all_DE_gene_annotation["gene_annotation"] <- diff_res_flu_all_annots_to_add
# write.csv(res_zaas_all_DE_gene_annotation, file = "../../../data/RNA-seq_data_chip/original_chip_data/DE_comparisons/res_zaas_all_annotations.csv")

# zaas specific 
res_zaas_specific_DE_gene_annotation <- diff_res_flu
diff_res_flu_annots_to_add <- gkey_all$GENENAME[match(diff_res_flu$gene_symbol, gkey_all$SYMBOL)]
res_zaas_specific_DE_gene_annotation["gene_annotation"] <- diff_res_flu_annots_to_add
# write.csv(res_zaas_specific_DE_gene_annotation, file = "../../../data/RNA-seq_data_chip/original_chip_data/DE_comparisons/res_zaas_specific_annotations.csv")
##--------------------------------------- Zaas et al GSEA (GAGE)
# enriched_paths_flu <- GSEA_Gage(xx1,DEmethod = "limma", G, time_point="","../../../figs/original_RNA-seq/pathway_analysis/GSEA_GAGE_Zaas_flu")
# enriched_paths_rsv <- GSEA_Gage(xx2,DEmethod = "limma", G, "../../../figs/original_RNA-seq/pathway_analysis/GSEA_GAGE_Zaas_rsv")
# enriched_paths_hsv <- GSEA_Gage(xx3,DEmethod = "limma", G, "../../../figs/original_RNA-seq/pathway_analysis/GSEA_GAGE_Zaas_hsv")
# separate viral groups ----
E_for_GAGE <- E
rownames(E_for_GAGE) <- G$ENTREZID[match(G$PROBEID, rownames(E))]
rhino_group <- which(S$virus.group == "Rhino")
rsv_group <- which(S$virus.group == "RSV")
flu_group <- which(S$virus.group == "H3N2")

symp <- which(S$clinical.group == 1)
asymp <- which(S$clinical.group == -1)

# tolerant are those that are asymptomatic but exhibit high viral shedding (high viral load)
# get samples at peak symptoms (if asymptomatic, jackson score < 6 over 5 days)
peak <- grep("peak",S$time.point)

flu_peak_symp <- intersect(flu_group,intersect(symp,peak))
flu_peak_asymp <- intersect(flu_group,intersect(asymp,peak))

rsv_peak_symp <- intersect(rsv_group,intersect(symp,peak))
rsv_peak_asymp <- intersect(rsv_group,intersect(asymp,peak))

hrv_peak_symp <- intersect(rhino_group,intersect(symp,peak))
hrv_peak_asymp <- intersect(rhino_group,intersect(asymp,peak))

enriched_paths_flu <- GSEA_Gage(E_for_GAGE, flu_peak_asymp, flu_peak_symp, gene_set_kegg = kegg.gs, output_loc="../../../figs/original_RNA-seq/pathway_analysis/GSEA_GAGE_Zaas_flu", time_point = "", num_paths = 12, draw_map = FALSE)
enriched_paths_rsv <- GSEA_Gage(E_for_GAGE, rsv_peak_asymp, rsv_peak_symp, gene_set_kegg = kegg.gs, output_loc="../../../figs/original_RNA-seq/pathway_analysis/GSEA_GAGE_Zaas_rsv", time_point = "", num_paths = 12, draw_map = FALSE)
enriched_paths_hrv <- GSEA_Gage(E_for_GAGE, hrv_peak_asymp, hrv_peak_symp, gene_set_kegg = kegg.gs, output_loc="../../../figs/original_RNA-seq/pathway_analysis/GSEA_GAGE_Zaas_hrv", time_point = "", num_paths = 12, draw_map = FALSE)


##--------------------------------------- Zaas et al GSEA (rpegeos) --- currently not working
flu_enr <-rpegeos_enrichment(diff_res_flu$gene_symbol, diff_res_flu, num_paths=20, method="limma",
                   org.Hs.eg.db = org.Hs.eg.db, output_loc="../../../figs/original_RNA-seq/pathway_analysis/GSEA_rpegeos_Zaas_flu")
## --------------------------------------

# get data
load("../../../data/Zaas_etal_microarray/2019-04-09_diffnet_zaas.RData")
load(file = "../../../data/2018-04-30_ttrust_data.RData")

network_data <- network_analysis(DN,ttrust_data, org.Hs.eg.db, gene_ids_chip, virus_specific, lung_genes)
node_names <- network_data$node_names
tmp <- network_data$tmp

# write files to load into cytoscape ----
# maybe do something for igraph instead? 
# write_delim(x = tmp, "res/2019-04-08_diffnet_results.txt", delim = "\t", col_names = TRUE)
# write_delim(x = node_names, "res/2019-04-08_diffnet_results_nodenames.txt", delim = "\t", col_names = TRUE)


# enrichment of nodes in subnetwork ----
subnet_gset <- cbind(node_names$ENTREZID, node_names$log2fc, node_names$p_val)
subnet_gset <- apply(subnet_gset, 2, as.numeric)
subnet_enr <- rpegeos::enrich_geneset(gene_set = subnet_gset)

subnet_enr %>% dplyr::filter(., probability_random < 0.01)

tmp %>% 
  #dplyr::filter(., x_lung == 1 | y_lung == 1) %>% 
  dplyr::filter(., x_tf == 1 | y_tf == 1) %>%
  dplyr::select(., x_name, y_name, x_tf, y_tf, cor1, cor2, change_type)

#--------------------------------------- comparing pathways (GAGE)
output_dir_compare_zaas_chip = "../../../figs/original_RNA-seq/pathway_analysis/"
out_name_compare_18hrs <-   "GSEA_GAGE_comparison_chip_zaas_18hrs"
out_name_compare_48hrs <- "GSEA_GAGE_comparison_chip_zaas_48hrs"
compare_chip_Zaas(enriched_paths_chip_18hrs, enriched_paths_flu, direction="upregulated", num_paths=15,out_name_compare_18hrs)
compare_chip_Zaas(enriched_paths_chip_48hrs, enriched_paths_flu, direction="upregulated", num_paths=15,out_name_compare_48hrs)

#--------------------------------------- comparing pathways (repgeos)
output_dir_compare_rpegeos_zaas_chip = "../../../figs/original_RNA-seq/pathway_analysis/"
compare_zaas_chip_rpegeos(flu_enr, chip_enr, numPaths = 15, output_dir_compare_rpegeos_zaas_chip)

## ---------------------------- Lung-specific genes
# plot ----
b2_18 <- res_v18$log2FoldChange[which(gene_ids_chip %in% lung_flu_biomarkers_18hrs)] # <-- chip data
#b3 <- y1$logFC[which(gene_mappings_copd$SYMBOL %in% b1)] # <-- copd data
b4_18 <- xx1$logFC[which(G$SYMBOL %in% lung_flu_biomarkers_18hrs)] # <-- zaas data, flu
b5_18 <- xx2$logFC[which(G$SYMBOL %in% lung_flu_biomarkers_18hrs)] # <-- zaas data, rsv
b6_18 <- xx3$logFC[which(G$SYMBOL %in% lung_flu_biomarkers_18hrs)] # <-- zaas data, hrv

dge_groups <- tibble(group = c(rep("this study", length(b2_18)),
                               #rep("copd data (benam et al)", length(b3)),
                               rep("flu infection data (zaas et al)", length(b4_18)),
                               rep("rsv infection data (zaas et al)", length(b5_18)),
                               rep("hrv infection data (zaas et al)", length(b6_18))),
                     gene_name = c(gene_ids_chip[gene_ids_chip %in% lung_flu_biomarkers_18hrs], 
                                   #gene_mappings_copd$SYMBOL[which(gene_mappings_copd$SYMBOL %in% b1)],
                                   G$SYMBOL[G$SYMBOL %in% lung_flu_biomarkers_18hrs],
                                   G$SYMBOL[G$SYMBOL %in% lung_flu_biomarkers_18hrs],
                                   G$SYMBOL[G$SYMBOL %in% lung_flu_biomarkers_18hrs]),
                     fold_data = c(b2_18, b4_18, b5_18, b6_18))

dge_groups %>% 
  dplyr::filter(., group != "copd data (benam et al)", fold_data > 0.5) %>%
  ggplot() + 
  geom_point(aes(x = group, y = fold_data, color = group), size = 5, alpha = 0.5) + 
  # scale_color_manual(values = c("purple", "orange")) +
  facet_wrap(. ~ gene_name) + 
  labs(x = NULL, y = "log2(fold difference to control)") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 12, color = "black", angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        legend.position = "none")

ggsave(filename = "../../../figs/original_RNA-seq/tissue_specificity/lung_specific_biomarkers_18hrs.pdf",
       plot = last_plot(), 
       device = "pdf", 
       height = 8,
       width = 10,
       scale = 1, 
       dpi = 400)

ggsave(filename = "../../../figs/original_RNA-seq/tissue_specificity/lung_specific_biomarkers_18hrs.png",
       plot = last_plot(), 
       device = "png", 
       height = 8,
       width = 10,
       scale = 1, 
       dpi = 400)

b2_48 <- res_v18$log2FoldChange[which(gene_ids_chip %in% lung_flu_biomarkers_48hrs)] # <-- chip data
#b3 <- y1$logFC[which(gene_mappings_copd$SYMBOL %in% b1)] # <-- copd data
b4_48 <- xx1$logFC[which(G$SYMBOL %in% lung_flu_biomarkers_48hrs)] # <-- zaas data, flu
b5_48 <- xx2$logFC[which(G$SYMBOL %in% lung_flu_biomarkers_48hrs)] # <-- zaas data, rsv
b6_48 <- xx3$logFC[which(G$SYMBOL %in% lung_flu_biomarkers_48hrs)] # <-- zaas data, hrv

dge_groups <- tibble(group = c(rep("this study", length(b2_48)),
                               #rep("copd data (benam et al)", length(b3)),
                               rep("flu infection data (zaas et al)", length(b4_48)),
                               rep("rsv infection data (zaas et al)", length(b5_48)),
                               rep("hrv infection data (zaas et al)", length(b6_48))),
                     gene_name = c(gene_ids_chip[gene_ids_chip %in% lung_flu_biomarkers_48hrs], 
                                   #gene_mappings_copd$SYMBOL[which(gene_mappings_copd$SYMBOL %in% b1)],
                                   G$SYMBOL[G$SYMBOL %in% lung_flu_biomarkers_48hrs],
                                   G$SYMBOL[G$SYMBOL %in% lung_flu_biomarkers_48hrs],
                                   G$SYMBOL[G$SYMBOL %in% lung_flu_biomarkers_48hrs]),
                     fold_data = c(b2_48, b4_48, b5_48, b6_48))

dge_groups %>% 
  dplyr::filter(., group != "copd data (benam et al)", fold_data > 0.5) %>%
  ggplot() + 
  geom_point(aes(x = group, y = fold_data, color = group), size = 5, alpha = 0.5) + 
  # scale_color_manual(values = c("purple", "orange")) +
  facet_wrap(. ~ gene_name) + 
  labs(x = NULL, y = "log2(fold difference to control)") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 12, color = "black", angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        legend.position = "none")

ggsave(filename = "../../../figs/original_RNA-seq/tissue_specificity/lung_specific_biomarkers_48hrs.pdf",
       plot = last_plot(), 
       device = "pdf", 
       height = 8,
       width = 10,
       scale = 1, 
       dpi = 400)

ggsave(filename = "../../../figs/original_RNA-seq/tissue_specificity/lung_specific_biomarkers_48hrs.png",
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

# differentially abundant signal peptides in our data
virus_sig_peps_18 <- intersect(sig_peps_names$SYMBOL, virus_specific_genes_18hrs)
virus_sig_peps_48 <- intersect(sig_peps_names$SYMBOL, virus_specific_genes_48hrs)

# secreted proteins ----
# based on HPA predicted data
hpa_pred <- read_delim("../../../data/secretome/2019-03-07_HPA_predicted_secreted_proteins.tsv", delim = "\t")
hpa_sec <- hpa_pred[grep("secreted", hpa_pred$`Protein class`), ]

# which virus secreted proteins are differentially abundant 
virus_secp_18 <- intersect(hpa_pred$Gene, virus_specific_genes_18hrs)
virus_secp_48 <- intersect(hpa_pred$Gene, virus_specific_genes_48hrs)
# predicted secreted biomarkers ----
virus_biom_sec_18 <- union(virus_secp_18, virus_sig_peps_18)
length(virus_biom_sec_18)
# virus_biom_sec # proteins that are predicted to be secreted and have signal peptides

# write.table(virus_biom_sec, file="../data/gene_lists/virus_secreted_signal_chip.txt", sep="\t", col.names = F, row.names = F, quote=F)
# write.table(virus_secp_48hrs, file="../../data/gene_lists/virus_secreted_signal_chip_48hrs.txt", sep="\t", col.names = F, row.names = F, quote=F)

virus_biom_sec_48 <- union(virus_secp_48, virus_sig_peps_48)

# secreted virus specific ----
d2_18 <- res_v18$log2FoldChange[which(gene_ids_chip %in% virus_biom_sec_18)] # <-- chip data
d4_18<- xx1$logFC[which(G$SYMBOL %in% virus_biom_sec_18)] # <-- zaas data, flu
d5_18 <- xx2$logFC[which(G$SYMBOL %in% virus_biom_sec_18)] # <-- zaas data, rsv
d6_18 <- xx3$logFC[which(G$SYMBOL %in% virus_biom_sec_18)] # <-- zaas data, hrv

sec_genes_18 <- tibble(group = c(rep("this study", length(d2_18)),
                                 rep("flu infection data (zaas et al)", length(d4_18)),
                                 rep("rsv infection data (zaas et al)", length(d5_18)),
                                 rep("hrv infection data (zaas et al)", length(d6_18))),
                       gene_name = c(gene_ids_chip[which(gene_ids_chip %in% virus_biom_sec_18)], 
                                     G$SYMBOL[which(G$SYMBOL %in% virus_biom_sec_18)],
                                     G$SYMBOL[which(G$SYMBOL %in% virus_biom_sec_18)],
                                     G$SYMBOL[which(G$SYMBOL %in% virus_biom_sec_18)]),
                       fold_data = c(d2_18, d4_18, d5_18, d6_18))

sec_genes_18 %>% 
  ggplot() + 
  geom_point(aes(x = forcats::fct_reorder(gene_name, fold_data, .desc = TRUE), y = fold_data, color = group), size = 3, alpha = 0.5) + 
  geom_hline(yintercept = 1, color = "black", lty = 2) + 
  # scale_color_manual(values = c("purple", "orange")) +
  labs(x = "Secreted biomarker", y = "log2(fold difference to control)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "bottom")

ggsave(filename = "../../../figs/original_RNA-seq/differential_expression/secreted_biomarkers_18hrs.pdf",
       plot = last_plot(), 
       device = "pdf", 
       height = 8,
       width = 10,
       scale = 1, 
       dpi = 400)

ggsave(filename = "../../../figs/original_RNA-seq/differential_expression/secreted_biomarkers_18hrs.png",
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
  geom_hline(yintercept = 1, color = "black", lty = 2) + 
  # scale_color_manual(values = c("purple", "orange")) +
  labs(x = "Secreted biomarker", y = "log2(fold difference to control)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
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

##--------------------------------------- Overlap analysis of feature importances
lin_svm_feats <- data.frame(data.table::fread('../data/models_analysis/SVM_linear_chip_SS_all_weights_annotation.csv', header=TRUE, sep = ','), row.names=1)
RF_feats <- data.frame(data.table::fread('../data/models_analysis/RF_chip_SS_training_all_weights_annotation.csv', header=TRUE, sep = ','), row.names=1)

lin_svm_feats = data.frame(rownames(lin_svm_feats), lin_svm_feats$Feature.importance)
RF_feats = data.frame(rownames(RF_feats), RF_feats$Feature.importance)


RRHO.ML <- RRHO(lin_svm_feats, RF_feats, stepsize = 2, alternative='enrichment', labels=c('linear SVM', 'RandomForest'), plots=TRUE, outputdir = '../figs')
pval <- pvalRRHO(RRHO.ML, 100)

