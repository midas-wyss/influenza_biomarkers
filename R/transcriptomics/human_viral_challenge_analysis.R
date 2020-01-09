##--------------------------------------- Influenza biomarkers project overview
# The goal of this project is to computationally identify biomkarers that validate the efficacy of an in-vitro "organ-on-chip" model of influenza infection. 
# Obtaining these biomakers will aid in i) evaluating the model's ability to recapitulate disease physiology ii) identifying novel clinical biomarkers that 
# distinguish influenza from other respiratory diseases and iii) aid in the development of therapeutics tageting influenza. These biomarkers will be recovered 
# from RNA-sequencing data capturing gene expression dynamics in a healthy and infected (COPD) model.  This project is being conducted in collaboration with the 
# Ingber lab and will utilize their acquired RNA-seq data alongside publicly available datasets. 

# Analysis of external data sets

# Author: Miguel A. Alcantar 
# Last updated: 07/12/2019

##--------------------------------------- Loading libraries needed for analysis

# If starting from scratch, you will likely need to install the following packages using BiocManager. 
#if (!requireNamespace("BiocManager")) --- likely won't need to run this
#install.packages("BiocManager")
#BiocManager::install("DESeq2")
#BiocManager::install("apeglm")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DiffNet")
#BiocManager::install("hgu133plus2.db")
#BiocManager::install("limma")
#BiocManager::install("GEOquery")
#BiocManager::install("affy")
#BiocManager::install("ArrayExpress")

# loading packages
library(DESeq2) # differential expression 
library(tidyverse) # for organizing data structures
library(readr) # for data parsing / loading in data files
library(readxl) # for reading Excel files
library(org.Hs.eg.db) # gene database
library(rpegeos) # created by Diogo (https://rdrr.io/github/diogocamacho/rpegeos/)
library(DiffNet)
library(hgu133plus2.db, quietly = TRUE) # human genome
library(limma, quietly = TRUE) # another differential expression package
library(affy)
library(oligo)
# packages for GSEA
library(pathview)
library(gage)
library(gageData)
library(ArrayExpress)
data(kegg.sets.hs) # KEGG pathways
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]

#library(DRUID) created by Diogo -- don't need yet (https://github.com/diogocamacho/druid/)
#library(grove)


source("load_RNAseq_transcriptomics_data.R") # script for loading data
source("DE_pathway_functions.R") # script for gene- and pathway-level analysis
source("preprocess_RNAseq_transcriptomics_data.R") # scripts for preprocessing data



##--------------------------------------- Huang et al analysis
load("../data/External_Human_Challenge_transcriptomics/External_microarray_Rdata/huang_etal_microarray.Rdata")

huang_E <- huang_data_all$expression
huang_G <- huang_data_all$genes
huang_S <- huang_data_all$metadata

# switching labels will just give you everything reversed
huang_DE_results <- DE_Huang_Woods(huang_S, huang_E)
GSEA_Gage(huang_DE_results,DEmethod = "limma", huang_G, "../figs/GSEA_GAGE_huang_etal")


diff_genes_huang <- which(huang_DE_results$adj.P.Val < 1 & abs(huang_DE_results$logFC) > 1)

# flu biomarkers ----
# these are genes that are differentially expressed in flu but not
# in the other 2 infection models

diff_res_huang <- tibble::tibble(gene_id=huang_G$ENTREZID[diff_genes_huang],
                               gene_symbol=huang_G$SYMBOL[diff_genes_huang],
                               fold_change=huang_DE_results$logFC[diff_genes_huang],
                               fdr_pvalue=huang_DE_results$adj.P.Val[diff_genes_huang])

huang_enr <-rpegeos_enrichment(diff_res_huang$gene_symbol, diff_res_huang, num_paths=20, method="limma",
                             org.Hs.eg.db = org.Hs.eg.db, output_loc="../figs/GSEA_rpegeos_huang")

### Temporal analysis

for(time_point in unique(huang_S$time_point)){
  idx_of_interst <- which(huang_S$time_point == time_point)
  huang_S_temp <- huang_S[idx_of_interst,]
  huang_E_temp <- huang_E[,idx_of_interst]
  
  huang_DE_results_temp <- DE_Huang_Woods(huang_S_temp, huang_E_temp)
  GSEA_Gage(huang_DE_results_temp, DEmethod = "limma", huang_G, paste("../figs/GSEA_GAGE_huang_etal", time_point, sep = "_"))
  
}



##--------------------------------------- Woods et al analysis
load("../data/External_Human_Challenge_transcriptomics/External_microarray_Rdata/woods_etal_microarray.Rdata")

woods_E <- woods_data_all$expression
woods_G <- woods_data_all$genes
woods_S <- woods_data_all$metadata

woods_DE_results <-DE_Huang_Woods(woods_S, woods_E)
GSEA_Gage(woods_DE_results,DEmethod = "limma", woods_G, "../figs/GSEA_GAGE_woods_etal")

diff_genes_woods <- which(woods_DE_results$adj.P.Val < 0.01 & abs(woods_DE_results$logFC) > 1)

# flu biomarkers ----
# these are genes that are differentially expressed in flu but not
# in the other 2 infection models

diff_res_woods <- tibble::tibble(gene_id=woods_G$ENTREZID[diff_genes_woods],
                                 gene_symbol=woods_G$SYMBOL[diff_genes_woods],
                                 fold_change=woods_DE_results$logFC[diff_genes_woods],
                                 fdr_pvalue=woods_DE_results$adj.P.Val[diff_genes_woods])

woods_enr <-rpegeos_enrichment(diff_res_woods$gene_symbol, diff_res_woods, num_paths=20, method="limma",
                               org.Hs.eg.db = org.Hs.eg.db, output_loc="../figs/GSEA_rpegeos_woods")
