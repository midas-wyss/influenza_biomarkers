##--------------------------------------- Influenza biomarkers project overview
# The goal of this project is to computationally identify biomkarers that validate the efficacy of an in-vitro "organ-on-chip" model of influenza infection. 
# Obtaining these biomakers will aid in i) evaluating the model's ability to recapitulate disease physiology ii) identifying novel clinical biomarkers that 
# distinguish influenza from other respiratory diseases and iii) aid in the development of therapeutics tageting influenza. These biomarkers will be recovered 
# from RNA-sequencing data capturing gene expression dynamics in a healthy and infected (COPD) model.  This project is being conducted in collaboration with the 
# Ingber lab and will utilize their acquired RNA-seq data alongside publicly available datasets. 

# RNA-seq is of lung epithelial cells

# This script will be used to conduct analyses on external Huang et al. microarray data
# taken from human viral challenge studies.

# Author: Miguel A. Alcantar 
# Last updated: 07/12/2019

##--------------------------------------- Installing packages
#BiocManager::install("pd.hugene.1.0.st.v1")
#BiocManager::install("hugene10sttranscriptcluster.db")
#BiocManager::install("hgu133a2.db")

# loading packages
library(tidyverse) # for organizing data structures
library(readr) # for data parsing / loading in data files
library(readxl) # for reading Excel files
library(org.Hs.eg.db) # gene database
library(DiffNet)
library(hgu133plus2.db, quietly = TRUE) # human genome
library(limma, quietly = TRUE) # another differential expression package
library(affy)
library(oligo)
library(pheatmap)
# packages for GSEA
library(pathview)
library(gage)
library(gageData)
library(ArrayExpress)
data(kegg.sets.hs) # KEGG pathways
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]

library("hugene10sttranscriptcluster.db")
library(hgu133a2.db) #  hgu133a2
library(illuminaHumanv4.db)

library("pd.hugene.1.0.st.v1")

##--------------------------------------- pre-processing Huang et al microarray data 
# (https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-30550/?query=Yongsheng+Huang&organism=Homo+sapiens&exptype%5B%5D=%22rna+assay%22)

# creating temporary directory that will contain microarray data
huang_etal_data_dir <- tempdir() 

if (!dir.exists(huang_etal_data_dir)) {
  dir.create(huang_etal_data_dir)
}

# downloading microaray data from ArrayExpress
anno_AE <- getAE("E-GEOD-30550", path = huang_etal_data_dir, type = "raw")

# reading in annotated data from Huang et al. -- SDRF gives mapping from sample to data
huang_sdrf_location <- file.path(huang_etal_data_dir, "E-GEOD-30550.sdrf.txt")
huang_SDRF <- read.delim(huang_sdrf_location)
rownames(huang_SDRF) <- huang_SDRF$Array.Data.File
huang_SDRF <- AnnotatedDataFrame(huang_SDRF)
huang_raw_data <- oligo::read.celfiles(filenames = file.path(huang_etal_data_dir,
                                                             huang_SDRF$Array.Data.File),
                                       verbose = FALSE, phenoData = huang_SDRF)
stopifnot(validObject(huang_raw_data))

# extracting phenotype data -- essentially metadata we care about
Biobase::pData(huang_raw_data) <- Biobase::pData(huang_raw_data)[, c("Source.Name", # GSM ID
                                                                     "Comment..Sample_description.", # RNA source + subject + time
                                                                     "Comment..Sample_source_name.", # RNA source + subject + time
                                                                     "Characteristics.clinic_pheno.", # A vs S
                                                                     "Characteristics.time_hpi.", # time point
                                                                     "FactorValue..CLINIC_PHENO.")] # A vs. S
# loading actual expression data
exp_raw <- log2(Biobase::exprs(huang_raw_data))
# PCA_raw <- prcomp(t(exp_raw), scale. = FALSE) -- this PCA takes a while

# boxplot to see distribution of expression for all samples -- there appear to be two outliers (consider removing)
# oligo::boxplot(huang_raw_data, target = "core",
#                main = "Boxplot of log2-intensitites for the raw data")

# subracting transcript median intensities to facilitate finding outliers
huang_eset <- oligo::rma(huang_raw_data, normalize = FALSE)
row_medians_assayData <- 
  Biobase::rowMedians(as.matrix(Biobase::exprs(huang_eset)))
RLE_data <- sweep(Biobase::exprs(huang_eset), 1, row_medians_assayData)
RLE_data <- as.data.frame(RLE_data)

RLE_data_gathered <-
  tidyr::gather(RLE_data, patient_array, log2_expression_deviation)

# boxplot of meadian-substracted samples -- still seem to be a few outliers
# ggplot2::ggplot(RLE_data_gathered, aes(patient_array,
#                                        log2_expression_deviation)) + 
#   geom_boxplot(outlier.shape = NA) + 
#   ylim(c(-2, 2)) + 
#   theme(axis.text.x = element_text(colour = "aquamarine4",
#                                    angle = 60, size = 6.5, hjust = 1 ,
#                                    face = "bold"))

# normalizing data
huang_eset_norm <- oligo::rma(huang_raw_data)

exp_huang <- Biobase::exprs(huang_eset_norm)

# pca on normalized data -- two outliers
PCA <- prcomp(t(exp_huang), scale = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     Disease =
                       Biobase::pData(huang_eset_norm)$FactorValue..CLINIC_PHENO.)

# ggplot(dataGG, aes(PC1, PC2)) +
#   geom_point(aes(shape = Disease)) +
#   ggtitle("PCA plot of the calibrated, summarized data") +
#   xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
#   ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   coord_fixed(ratio = sd_ratio) +
#   scale_shape_manual(values = c(4,15)) +
#   scale_color_manual(values = c("darkorange2", "dodgerblue4"))

disease_names <- ifelse(str_detect(pData
                                   (huang_eset_norm)$FactorValue..CLINIC_PHENO.,
                                   "Asympto"), "A", "S")
# heatmap of data -- there does seem to be a little bit of clustering for Asymptomatic and Symptomatic patients
# consider checking if they cluster based on severity of symptoms  
annotation_for_heatmap <-
  data.frame(Disease = disease_names)

row.names(annotation_for_heatmap) <- row.names(pData(huang_eset_norm))

# dists <- as.matrix(dist(t(exp_huang), method = "euclidean"))
# 
# rownames(dists) <- row.names(pData(huang_eset_norm))
# hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(255))
# colnames(dists) <- NULL
# diag(dists) <- NA
# 
# ann_colors <- list(Disease = c(A = "blue4", S = "cadetblue2"))
# 
# pheatmap(dists, col = (hmcol),
#          annotation_row = annotation_for_heatmap,
#          annotation_colors = ann_colors,
#          legend = TRUE,
#          treeheight_row = 0,
#          legend_breaks = c(min(dists, na.rm = TRUE),
#                            max(dists, na.rm = TRUE)),
#          legend_labels = (c("small distance", "large distance")),
#          main = "Clustering heatmap for the calibrated samples")                   
huang_medians <- rowMedians(Biobase::exprs(huang_eset_norm))

# distribution of microarray intensities
hist_res <- hist(huang_medians, 100, col = "cornsilk1", freq = FALSE,
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")
man_threshold <- 4 # intensity threshold

hist_res <- hist(huang_medians, 100, col = "cornsilk", freq = FALSE,
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")

abline(v = man_threshold, col = "coral4", lwd = 2)


# removing probes that fall below the threshold for #samples > size if smalles group
no_of_samples <-
  table(paste0(pData(huang_eset_norm)$FactorValue..CLINIC_PHENO.))
samples_cutoff <- min(no_of_samples)

idx_man_threshold <- apply(Biobase::exprs(huang_eset_norm), 1,
                           function(x){
                             sum(x > man_threshold) >= samples_cutoff})
table(idx_man_threshold)

huang_manfiltered <- subset(huang_eset_norm, idx_man_threshold)

# mapping probe ID to gene name -- the hgu133a2.db depends on the platform used
anno_huang <-AnnotationDbi::select(hgu133a2.db,
                      keys = featureNames(huang_manfiltered),
                      columns = c("SYMBOL", "GENENAME"),
                      keytype = "PROBEID")

anno_huang <- subset(anno_huang, !is.na(SYMBOL))

anno_grouped <- group_by(anno_huang, PROBEID)
anno_summarized <-
  dplyr::summarize(anno_grouped, no_of_matches = n_distinct(SYMBOL))
head(anno_summarized)

# filtering data
anno_filtered <- filter(anno_summarized, no_of_matches > 1)

probe_stats <- anno_filtered
ids_to_exlude <- (featureNames(huang_manfiltered) %in% probe_stats$PROBEID)
# table(ids_to_exlude)
huang_final <- subset(huang_manfiltered, !ids_to_exlude)
fData(huang_final)$PROBEID <- rownames(fData(huang_final))
rownames(fData(huang_final)) <- fData(huang_final)$PROBEID
# validObject(palmieri_final)

# saving expression data before subsequent filtering
Huang_etal_expression_data <- exprs(huang_final)

# getting gene data only for filtered genes
anno_huang_filtered <-AnnotationDbi::select(hgu133a2.db,
                                      keys = featureNames(huang_final),
                                      columns = c("SYMBOL", "GENENAME"),
                                      keytype = "PROBEID")
geneIDs_huang <- anno_huang_filtered$SYMBOL

# mapping to ENTREZ IDs
gkey_huang <- AnnotationDbi::select(x = org.Hs.eg.db, 
                              keys = geneIDs_huang, 
                              columns = "ENTREZID", 
                              keytype = "SYMBOL")

# removing genes with multiple ENTREZ IDs
nix <- which(sapply(anno_huang_filtered$SYMBOL, function(y) length(which(gkey_huang$SYMBOL == y))) != 1)
if(length(nix) != 0) {
  geneIDs_huang <- geneIDs_huang[-nix]
  Huang_etal_expression_data <- Huang_etal_expression_data[-nix, ]
}

# final expression matrix
Huang_etal_expression_data_final <- Huang_etal_expression_data

# grabbing genes with only 1 ENTREZ ID mapping
gene_data_pre <- subset(anno_huang_filtered, anno_huang_filtered$PROBEID %in% rownames(Huang_etal_expression_data))

gkey_huangv3 <- AnnotationDbi::select(x = org.Hs.eg.db, 
                                      keys = gene_data_pre$SYMBOL, 
                                      columns = "ENTREZID", 
                                      keytype = "SYMBOL")
# final genes after preprocessing
gene_data_final <- gene_data_pre %>% 
  tibble::add_column(., ENTREZID = gkey_huangv3$ENTREZID, .before = 1)

huang_subject_IDs <- vector(mode="character", length=length(levels(huang_final$Comment..Sample_source_name)))
huang_time_point <- vector(mode="character", length=length(levels(huang_final$Comment..Sample_source_name)))
counter_huang = 1

# finding time and corresponding sample ID
for (value in levels(huang_final$Comment..Sample_source_name)){
  huang_subject_IDs[[counter_huang]] <-  substr(unlist(strsplit(value[[1]], split = " "))[[5]],1,2)
  if(length(unlist(strsplit(value[[1]], split = " "))) == 6 )
  {
    huang_time_point[[counter_huang]] <-  unlist(strsplit(value[[1]], split = " "))[[6]]
  }
  else{
    huang_time_point[[counter_huang]] <-  unlist(strsplit(value[[1]], split = " "))[[7]]
  }
  
  counter_huang = counter_huang + 1
  # print(value)
}

# final metadata
huang_metadata <- tibble::tibble(source_name=huang_final$Source.Name,
                               disease=huang_final$FactorValue..CLINIC_PHENO.,
                               time_point=huang_time_point,
                               subject_ID = huang_subject_IDs)
# list with all the important Huang et al. data
huang_data_all <- list("expression" = Huang_etal_expression_data_final, "metadata" = huang_metadata, "genes" = gene_data_final)

save(huang_data_all, file = "../data/External_Human_Challenge_transcriptomics/External_microarray_Rdata/huang_etal_microarray.Rdata")


### --------------------------------------- woods et al micro array
# https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-52428/?query=Geoffrey+S+Ginsburg&organism=Homo+sapiens&exptype%5B%5D=%22rna+assay%22

# creating temporary directory that will contain microarray data from woods et al.
woods_etal_data_dir <- tempdir() 

if (!dir.exists(woods_etal_data_dir)) {
  dir.create(woods_etal_data_dir)
}

# downloading microaray data from ArrayExpress
woods_anno_AE <- getAE("E-GEOD-52428", path = woods_etal_data_dir, type = "raw")

# reading in annotated data from Huang et al. -- SDRF gives mapping from sample to data
woods_sdrf_location <- file.path(woods_etal_data_dir, "E-GEOD-52428.sdrf.txt")
woods_SDRF <- read.delim(woods_sdrf_location)
rownames(woods_SDRF) <- woods_SDRF$Array.Data.File
woods_SDRF <- AnnotatedDataFrame(woods_SDRF)
woods_raw_data <- oligo::read.celfiles(filenames = file.path(woods_etal_data_dir,
                                                             woods_SDRF$Array.Data.File),
                                       verbose = FALSE, phenoData = woods_SDRF)
stopifnot(validObject(woods_raw_data))



# extracting phenotype data -- essentially metadata we care about
Biobase::pData(woods_raw_data) <- Biobase::pData(woods_raw_data)[, c("Source.Name", # GSM ID
                                                                     "Characteristics..challenge.virus.", # RNA source + subject + time
                                                                     "Characteristics..clinical.phenotype.", # RNA source + subject + time
                                                                     "Characteristics..subject.id.", # A vs S
                                                                     "Characteristics..total.symptom.score.", # time point
                                                                     "FactorValue..TIME.")] # A vs. S
# loading actual expression data
woods_exp_raw <- log2(Biobase::exprs(woods_raw_data))
# PCA_raw <- prcomp(t(exp_raw), scale. = FALSE) -- this PCA takes a while

# boxplot to see distribution of expression for all samples -- there appear to be two outliers (consider removing)
# oligo::boxplot(woods_raw_data, target = "core",
              #  main = "Boxplot of log2-intensitites for the raw data")

# subracting transcript median intensities to facilitate finding outliers
woods_eset <- oligo::rma(woods_raw_data, normalize = FALSE)
woods_row_medians_assayData <- 
  Biobase::rowMedians(as.matrix(Biobase::exprs(woods_eset)))
woods_RLE_data <- sweep(Biobase::exprs(woods_eset), 1, woods_row_medians_assayData)
woods_RLE_data <- as.data.frame(woods_RLE_data)

woods_RLE_data_gathered <-
  tidyr::gather(woods_RLE_data, patient_array, log2_expression_deviation)

# boxplot of meadian-substracted samples -- still seem to be a few outliers
# ggplot2::ggplot(woods_RLE_data_gathered, aes(patient_array,
  #                                      log2_expression_deviation)) + 
  # geom_boxplot(outlier.shape = NA) + 
  # ylim(c(-2, 2)) + 
  # theme(axis.text.x = element_text(colour = "aquamarine4",
  #                                  angle = 60, size = 6.5, hjust = 1 ,
  #                                  face = "bold"))

# normalizing data
woods_eset_norm <- oligo::rma(woods_raw_data)

exp_woods <- Biobase::exprs(woods_eset_norm)

# pca on normalized data -- two outliers
PCA <- prcomp(t(exp_woods), scale = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

woods_dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     Disease =
                       Biobase::pData(woods_eset_norm)$Characteristics..clinical.phenotype.)

ggplot(woods_dataGG, aes(PC1, PC2)) +
  geom_point(aes(shape = Disease)) +
  ggtitle("PCA plot of the calibrated, summarized data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(4,15)) +
  scale_color_manual(values = c("darkorange2", "dodgerblue4"))

woods_disease_names <- ifelse(str_detect(pData
                                   (woods_eset_norm)$Characteristics..clinical.phenotype.,
                                   "Asympto"), "A", "S")
# heatmap of data -- there does seem to be a little bit of clustering for Asymptomatic and Symptomatic patients
# consider checking if they cluster based on severity of symptoms  
woods_annotation_for_heatmap <-
  data.frame(Disease = woods_disease_names)

row.names(woods_annotation_for_heatmap) <- row.names(pData(woods_eset_norm))

dists_woods <- as.matrix(dist(t(exp_woods), method = "euclidean"))

rownames(dists_woods) <- row.names(pData(woods_eset_norm))
hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(255))
colnames(dists_woods) <- NULL
diag(dists_woods) <- NA

ann_colors <- list(Disease = c(A = "blue4", S = "cadetblue2"))

pheatmap(dists_woods, col = (hmcol),
         annotation_row = woods_annotation_for_heatmap,
         annotation_colors = ann_colors,
         legend = TRUE,
         treeheight_row = 0,
         legend_breaks = c(min(dists, na.rm = TRUE),
                           max(dists, na.rm = TRUE)),
         legend_labels = (c("small distance", "large distance")),
         main = "Clustering heatmap for the calibrated samples")                   
woods_medians <- rowMedians(Biobase::exprs(woods_eset_norm))

# distribution of microarray intensities
woods_hist_res <- hist(woods_medians, 100, col = "cornsilk1", freq = FALSE,
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")
woods_man_threshold <- 4 # intensity threshold

woods_hist_res <- hist(woods_medians, 100, col = "cornsilk", freq = FALSE,
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")

abline(v = woods_man_threshold, col = "coral4", lwd = 2)


# removing probes that fall below the threshold for #samples > size if smalles group
woods_no_of_samples <-
  table(paste0(pData(woods_eset_norm)$Characteristics..clinical.phenotype.))
woods_samples_cutoff <- min(woods_no_of_samples)

woods_idx_man_threshold <- apply(Biobase::exprs(woods_eset_norm), 1,
                           function(x){
                             sum(x > woods_man_threshold) >= woods_samples_cutoff})
table(woods_idx_man_threshold)

woods_manfiltered <- subset(woods_eset_norm, woods_idx_man_threshold)

# mapping probe ID to gene name -- the illuminaHumanv4.db depends on the platform used
anno_woods <-AnnotationDbi::select(hgu133a2.db,
                                      keys = featureNames(woods_manfiltered),
                                      columns = c("SYMBOL", "GENENAME"),
                                      keytype = "PROBEID")

anno_woods <- subset(anno_woods, !is.na(SYMBOL))

anno_grouped_woods <- group_by(anno_woods, PROBEID)
anno_summarized_woods <-
  dplyr::summarize(anno_grouped_woods, no_of_matches = n_distinct(SYMBOL))
head(anno_summarized_woods)

# filtering data
anno_filtered_woods <- filter(anno_summarized_woods, no_of_matches > 1)

probe_stats_woods <- anno_filtered_woods
# head(probe_stats_woods)
# nrow(probe_stats_woods)
ids_to_exlude_woods <- (featureNames(woods_manfiltered) %in% probe_stats_woods$PROBEID)
# table(ids_to_exlude_woods)
woods_final <- subset(woods_manfiltered, !ids_to_exlude_woods)
# head(woods_final)
# validObject(woods_final)
#head(anno_woods)
fData(woods_final)$PROBEID <- rownames(fData(woods_final))
rownames(fData(woods_final)) <- fData(woods_final)$PROBEID
# validObject(woods_final)

# saving expression data before subsequent filtering
woods_etal_expression_data <- exprs(woods_final)

# getting gene data only for filtered genes
anno_woods_filtered <-AnnotationDbi::select(hgu133a2.db,
                                               keys = featureNames(woods_final),
                                               columns = c("SYMBOL", "GENENAME"),
                                               keytype = "PROBEID")
geneIDs_woods <- anno_woods_filtered$SYMBOL

# mapping to ENTREZ IDs
gkey_woods<- AnnotationDbi::select(x = org.Hs.eg.db, 
                                    keys = geneIDs_woods, 
                                    columns = "ENTREZID", 
                                    keytype = "SYMBOL")

# removing genes with multiple ENTREZ IDs
nix <- which(sapply(anno_woods_filtered$SYMBOL, function(y) length(which(gkey_woods$SYMBOL == y))) != 1)
if(length(nix) != 0) {
  geneIDs_woods <- geneIDs_woods[-nix]
  woods_etal_expression_data <- woods_etal_expression_data[-nix, ]
}

# final expression matrix
woods_etal_expression_data_final <- woods_etal_expression_data

# grabbing genes with only 1 ENTREZ ID mapping
gene_data_pre_woods <- subset(anno_woods_filtered, anno_woods_filtered$PROBEID %in% rownames(woods_etal_expression_data))

gkey_woodsv3 <- AnnotationDbi::select(x = org.Hs.eg.db, 
                                      keys = gene_data_pre_woods$SYMBOL, 
                                      columns = "ENTREZID", 
                                      keytype = "SYMBOL")
# final genes after preprocessing
gene_data_final_woods <- gene_data_pre_woods %>% 
  tibble::add_column(., ENTREZID = gkey_woodsv3$ENTREZID, .before = 1)

# final metadata
woods_metadata <- tibble::tibble(source_name=woods_final$Source.Name,
                                 disease=woods_final$Characteristics..clinical.phenotype.,
                                 time_point=woods_final$FactorValue..TIME.,
                                 virus = woods_final$Characteristics..challenge.virus.,
                                 subject_virus = woods_final$Characteristics..subject.id.)
# list with all the important Huang et al. data
woods_data_all <- list("expression" = woods_etal_expression_data_final, "metadata" = woods_metadata, "genes" = gene_data_final_woods)

save(woods_data_all, file = "../data/External_Human_Challenge_transcriptomics/External_microarray_Rdata/woods_etal_microarray.Rdata")
# PCA <- prcomp(t(woods_etal_expression_data), scale = FALSE)
# 
# percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
# sd_ratio <- sqrt(percentVar[2] / percentVar[1])
# 
# dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
#                      Disease =
#                        Biobase::pData(woods_final)$Characteristics..clinical.phenotype.)
# 
# ggplot(dataGG, aes(PC1, PC2)) +
#   geom_point(aes(color = Disease)) +
#   ggtitle("PCA plot of the calibrated, summarized data") +
#   xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
#   ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   coord_fixed(ratio = sd_ratio) +
#   scale_shape_manual(values = c(4,15)) +
#   scale_color_manual(values = c("darkorange2", "dodgerblue4"))
# 
# # boxplot to see distribution of expression for all samples -- there appear to be two outliers (consider removing)
# oligo::boxplot(woods_etal_expression_data, target = "core",
#                main = "Boxplot of log2-intensitites for the raw data")



### --------------------------------------- davenport 

# creating temporary directory that will contain microarray data from Davenport et al.
davenport_etal_data_dir <- tempdir() 

if (!dir.exists(davenport_etal_data_dir)) {
  dir.create(davenport_etal_data_dir)
}

# downloading microaray data from ArrayExpress
davenport_anno_AE <- getAE("E-GEOD-61754", path = davenport_etal_data_dir, type = "raw")

# reading in annotated data from davenport et al. -- SDRF gives mapping from sample to data
davenport_sdrf_location <- file.path(davenport_etal_data_dir, "E-GEOD-61754.sdrf.txt")
davenport_SDRF <- read.delim(davenport_sdrf_location)
rownames(davenport_SDRF) <- davenport_SDRF$Array.Data.File
davenport_SDRF <- AnnotatedDataFrame(davenport_SDRF)
davenport_raw_data <- oligo::read.celfiles(filenames = file.path(davenport_etal_data_dir,
                                                                 davenport_SDRF$Array.Data.File),
                                       verbose = FALSE, phenoData = davenport_SDRF)

stopifnot(validObject(davenport_SDRF))



# extracting phenotype data -- essentially metadata we care about
Biobase::pData(woods_raw_data) <- Biobase::pData(woods_raw_data)[, c("Source.Name", # GSM ID
                                                                     "Characteristics..challenge.virus.", # RNA source + subject + time
                                                                     "Characteristics..clinical.phenotype.", # RNA source + subject + time
                                                                     "Characteristics..subject.id.", # A vs S
                                                                     "Characteristics..total.symptom.score.", # time point
                                                                     "FactorValue..TIME.")] # A vs. S
# loading actual expression data
woods_exp_raw <- log2(Biobase::exprs(woods_raw_data))
# PCA_raw <- prcomp(t(exp_raw), scale. = FALSE) -- this PCA takes a while

# boxplot to see distribution of expression for all samples -- there appear to be two outliers (consider removing)
oligo::boxplot(woods_raw_data, target = "core",
               main = "Boxplot of log2-intensitites for the raw data")

# subracting transcript median intensities to facilitate finding outliers
woods_eset <- oligo::rma(woods_raw_data, normalize = FALSE)
woods_row_medians_assayData <- 
  Biobase::rowMedians(as.matrix(Biobase::exprs(woods_eset)))
woods_RLE_data <- sweep(Biobase::exprs(woods_eset), 1, woods_row_medians_assayData)
woods_RLE_data <- as.data.frame(woods_RLE_data)

woods_RLE_data_gathered <-
  tidyr::gather(woods_RLE_data, patient_array, log2_expression_deviation)

# boxplot of meadian-substracted samples -- still seem to be a few outliers
ggplot2::ggplot(woods_RLE_data_gathered, aes(patient_array,
                                             log2_expression_deviation)) + 
  geom_boxplot(outlier.shape = NA) + 
  ylim(c(-2, 2)) + 
  theme(axis.text.x = element_text(colour = "aquamarine4",
                                   angle = 60, size = 6.5, hjust = 1 ,
                                   face = "bold"))

# normalizing data
woods_eset_norm <- oligo::rma(woods_raw_data)

exp_woods <- Biobase::exprs(woods_eset_norm)

# pca on normalized data -- two outliers
PCA <- prcomp(t(exp_woods), scale = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

woods_dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                           Disease =
                             Biobase::pData(woods_eset_norm)$Characteristics..clinical.phenotype.)

ggplot(woods_dataGG, aes(PC1, PC2)) +
  geom_point(aes(shape = Disease)) +
  ggtitle("PCA plot of the calibrated, summarized data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(4,15)) +
  scale_color_manual(values = c("darkorange2", "dodgerblue4"))

woods_disease_names <- ifelse(str_detect(pData
                                         (woods_eset_norm)$Characteristics..clinical.phenotype.,
                                         "Asympto"), "A", "S")
# heatmap of data -- there does seem to be a little bit of clustering for Asymptomatic and Symptomatic patients
# consider checking if they cluster based on severity of symptoms  
woods_annotation_for_heatmap <-
  data.frame(Disease = woods_disease_names)

row.names(woods_annotation_for_heatmap) <- row.names(pData(woods_eset_norm))

dists_woods <- as.matrix(dist(t(exp_woods), method = "euclidean"))

rownames(dists_woods) <- row.names(pData(woods_eset_norm))
hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(255))
colnames(dists_woods) <- NULL
diag(dists_woods) <- NA

ann_colors <- list(Disease = c(A = "blue4", S = "cadetblue2"))

pheatmap(dists_woods, col = (hmcol),
         annotation_row = woods_annotation_for_heatmap,
         annotation_colors = ann_colors,
         legend = TRUE,
         treeheight_row = 0,
         legend_breaks = c(min(dists, na.rm = TRUE),
                           max(dists, na.rm = TRUE)),
         legend_labels = (c("small distance", "large distance")),
         main = "Clustering heatmap for the calibrated samples")                   
woods_medians <- rowMedians(Biobase::exprs(woods_eset_norm))

# distribution of microarray intensities
woods_hist_res <- hist(woods_medians, 100, col = "cornsilk1", freq = FALSE,
                       main = "Histogram of the median intensities",
                       border = "antiquewhite4",
                       xlab = "Median intensities")
woods_man_threshold <- 4 # intensity threshold

woods_hist_res <- hist(woods_medians, 100, col = "cornsilk", freq = FALSE,
                       main = "Histogram of the median intensities",
                       border = "antiquewhite4",
                       xlab = "Median intensities")

abline(v = woods_man_threshold, col = "coral4", lwd = 2)


# removing probes that fall below the threshold for #samples > size if smalles group
woods_no_of_samples <-
  table(paste0(pData(woods_eset_norm)$Characteristics..clinical.phenotype.))
woods_samples_cutoff <- min(woods_no_of_samples)

woods_idx_man_threshold <- apply(Biobase::exprs(woods_eset_norm), 1,
                                 function(x){
                                   sum(x > woods_man_threshold) >= woods_samples_cutoff})
table(woods_idx_man_threshold)

woods_manfiltered <- subset(woods_eset_norm, woods_idx_man_threshold)

# mapping probe ID to gene name -- the illuminaHumanv4.db depends on the platform used
anno_woods <-AnnotationDbi::select(hgu133a2.db,
                                   keys = featureNames(woods_manfiltered),
                                   columns = c("SYMBOL", "GENENAME"),
                                   keytype = "PROBEID")

anno_woods <- subset(anno_woods, !is.na(SYMBOL))

anno_grouped_woods <- group_by(anno_woods, PROBEID)
anno_summarized_woods <-
  dplyr::summarize(anno_grouped_woods, no_of_matches = n_distinct(SYMBOL))
head(anno_summarized_woods)

# filtering data
anno_filtered_woods <- filter(anno_summarized_woods, no_of_matches > 1)

probe_stats_woods <- anno_filtered_woods
# head(probe_stats_woods)
# nrow(probe_stats_woods)
ids_to_exlude_woods <- (featureNames(woods_manfiltered) %in% probe_stats_woods$PROBEID)
# table(ids_to_exlude_woods)
woods_final <- subset(woods_manfiltered, !ids_to_exlude_woods)
# head(woods_final)
# validObject(woods_final)
#head(anno_woods)
fData(woods_final)$PROBEID <- rownames(fData(woods_final))
rownames(fData(woods_final)) <- fData(woods_final)$PROBEID
# validObject(woods_final)

# saving expression data before subsequent filtering
woods_etal_expression_data <- exprs(woods_final)

# getting gene data only for filtered genes
anno_woods_filtered <-AnnotationDbi::select(hgu133a2.db,
                                            keys = featureNames(woods_final),
                                            columns = c("SYMBOL", "GENENAME"),
                                            keytype = "PROBEID")
geneIDs_woods <- anno_woods_filtered$SYMBOL

# mapping to ENTREZ IDs
gkey_woods<- AnnotationDbi::select(x = org.Hs.eg.db, 
                                   keys = geneIDs_woods, 
                                   columns = "ENTREZID", 
                                   keytype = "SYMBOL")

# removing genes with multiple ENTREZ IDs
nix <- which(sapply(anno_woods_filtered$SYMBOL, function(y) length(which(gkey_woods$SYMBOL == y))) != 1)
if(length(nix) != 0) {
  geneIDs_woods <- geneIDs_woods[-nix]
  woods_etal_expression_data <- woods_etal_expression_data[-nix, ]
}

# final expression matrix
woods_etal_expression_data_final <- woods_etal_expression_data

# grabbing genes with only 1 ENTREZ ID mapping
gene_data_pre_woods <- subset(anno_woods_filtered, anno_woods_filtered$PROBEID %in% rownames(woods_etal_expression_data))

gkey_woodsv3 <- AnnotationDbi::select(x = org.Hs.eg.db, 
                                      keys = gene_data_pre_woods$SYMBOL, 
                                      columns = "ENTREZID", 
                                      keytype = "SYMBOL")
# final genes after preprocessing
gene_data_final_woods <- gene_data_pre_woods %>% 
  tibble::add_column(., ENTREZID = gkey_woodsv3$ENTREZID, .before = 1)

# final metadata
woods_metadata <- tibble::tibble(source_name=woods_final$Source.Name,
                                 disease=woods_final$Characteristics..clinical.phenotype.,
                                 time_point=woods_final$FactorValue..TIME.)
# list with all the important Huang et al. data
woods_data_all <- list("expression" = davenport_etal_expression_data_final, "metadata" = davenport_metadata, "genes" = gene_data_final_davenport)

save(woods_data_all, file = "../data/External_Human_Challenge_transcriptomics/External_microarray_Rdata/woods_etal_microarray.Rdata")
# PCA <- prcomp(t(woods_etal_expression_data), scale = FALSE)
# 
# percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
# sd_ratio <- sqrt(percentVar[2] / percentVar[1])
# 
# dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
#                      Disease =
#                        Biobase::pData(woods_final)$Characteristics..clinical.phenotype.)
# 
# ggplot(dataGG, aes(PC1, PC2)) +
#   geom_point(aes(color = Disease)) +
#   ggtitle("PCA plot of the calibrated, summarized data") +
#   xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
#   ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   coord_fixed(ratio = sd_ratio) +
#   scale_shape_manual(values = c(4,15)) +
#   scale_color_manual(values = c("darkorange2", "dodgerblue4"))
# 
# # boxplot to see distribution of expression for all samples -- there appear to be two outliers (consider removing)
# oligo::boxplot(woods_etal_expression_data, target = "core",
#                main = "Boxplot of log2-intensitites for the raw data")



