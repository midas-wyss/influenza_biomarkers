# Transcriptomics work

# libraries
library(tidyverse)
library(readr)
library(readxl)
library(DRUID)
library(grove)
library(DESeq2)
library(org.Hs.eg.db)
library(rpegeos)
library(diffnet)

# functions
source("R/stouffer_specificity.R")
source("R/limma_dge.R")


# data from project ----
# count data 
count_data <- read_csv("data/rna_seq/2018-05-29_STAR_Gene_Counts.csv")
gene_ids <- count_data$Gene_ID
count_data <- count_data %>% dplyr::select(., -Gene_ID)

# metadata
metadata <- readxl::read_xlsx("data/rna_seq/2018-05-29_RNAseq_Sample IDs.xlsx")
t1 <- gsub("5/2", "May2", gsub("5/1", "May1", gsub("-", "_", metadata$`Tube ID`)))

metadata <- metadata %>% 
  dplyr::mutate(., sample_name = t1)

metadata <- metadata[match(colnames(count_data), metadata$sample_name), ]

group <- character(length = nrow(metadata))
group[metadata$Sample == "Control"] <- "control"
group[metadata$Sample == "Virus Treated"] <- "virus"
group[metadata$Sample == "Poly IC Treated"] <- "poly_ic"

metadata <- metadata %>% 
  dplyr::mutate(., group = group)

sample_data <- tibble(sample_name = metadata$sample_name,
                      time_point = as.numeric(gsub("h", "", metadata$`Time point`)),
                      group = metadata$group)

sample_data <- sample_data %>% 
  dplyr::mutate(., group_time = paste(group, time_point, sep = "_"))


# process data ahead of analysis ----
# remove probes with zero counts throughout
nix <- which(rowSums(count_data) == 0)
if(length(nix) != 0) {
  gene_ids <- gene_ids[-nix]
  count_data <- count_data[-nix, ]
}

# remove gene ids with multiple entrez mappings
gkey <- AnnotationDbi::select(x = org.Hs.eg.db, 
                              keys = gene_ids, 
                              columns = c("ENTREZID", "GENENAME"), 
                              keytype = "SYMBOL")

nix <- which(sapply(gene_ids, function(y) length(which(gkey$SYMBOL == y))) != 1)
if(length(nix) != 0) {
  gene_ids <- gene_ids[-nix]
  count_data <- count_data[-nix, ]
}

# remove genes with NA entrez id
gkey <- AnnotationDbi::select(x = org.Hs.eg.db, 
                              keys = gene_ids, 
                              columns = c("ENTREZID", "GENENAME"), 
                              keytype = "SYMBOL")

nix <- which(is.na(gkey$ENTREZID))
if(length(nix) != 0) {
  gene_ids <- gene_ids[-nix]
  count_data <- count_data[-nix, ]
}

# frequency across samples
# minimum 5 read counts, at least in 75% of the samples
x1 <- apply(count_data, 1, function(y) length(which(y >= 5)))
x1 <- which(x1 >= round(0.75 * ncol(count_data)))

count_data <- count_data[x1, ]
gene_ids <- gene_ids[x1]
gkey <- AnnotationDbi::select(x = org.Hs.eg.db, 
                              keys = gene_ids, 
                              columns = c("ENTREZID", "GENENAME"),
                              keytype = "SYMBOL")

# check for samples similarity
# a1 <- cor(count_data[, colnames(count_data) %in% sample_data$sample_name[which(sample_data$group_time == "virus_18")]])
# a2 <- cor(count_data[, colnames(count_data) %in% sample_data$sample_name[which(sample_data$group_time == "virus_48")]])
# a3 <- cor(count_data[, colnames(count_data) %in% sample_data$sample_name[which(sample_data$group_time == "control_18")]])
# a4 <- cor(count_data[, colnames(count_data) %in% sample_data$sample_name[which(sample_data$group_time == "control_48")]])
# a5 <- cor(count_data[, colnames(count_data) %in% sample_data$sample_name[which(sample_data$group_time == "poly_ic_18")]])

# samples to remove have correlation < 0.98 with all others
# ig_samples <- c("V1_May1", "V1_May2", "C1_May1") 
# 
# sample_data <- sample_data[-which(sample_data$sample_name %in% ig_samples), ]
# count_data <- count_data[, -which(colnames(count_data) %in% ig_samples)]

# differential expression ----
# use deseq2
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = sample_data,
                              design = ~ group_time)

dds <- DESeq(object = dds)

# ec_e_ab_res <- results(dds,contrast=c("Group","H_E_AB","M_E_AB"),alpha=0.05,lfcThreshold = 1)
contrast_v18 <- c("group_time", "virus_18", "control_18")
contrast_ic18 <- c("group_time", "poly_ic_18", "control_18")

res_v18 <- results(dds,
        contrast = contrast_v18,
        pAdjustMethod = "fdr",
        cooksCutoff = FALSE)

res_v18 <- res_v18 %>% 
  as_tibble() %>%
  tibble::add_column(., gene = gkey$SYMBOL, .before = 1) %>% 
  tibble::add_column(., entrez_id = gkey$ENTREZID, .before = 2) %>% 
  tibble::add_column(., sig = 0) %>%
  dplyr::mutate(., sig = replace(sig, list = which(padj < 0.01 & abs(log2FoldChange) > 1), 1))
  # dplyr::mutate(., sig = replace(sig, list = which(padj < 0.01), 1))

res_v18 %>% 
  ggplot() + 
  geom_point(aes(y = log2FoldChange, x = log10(baseMean), color = as.factor(sig)), alpha = 0.5, size = 3) + 
  ylim(c(-max(abs(res_v18$log2FoldChange)), max(abs(res_v18$log2FoldChange)))) +
  # scale_color_viridis_d() +
  scale_color_manual(breaks = c(0, 1), name = NULL, values = c("black", "red"), labels = c("Not significant", "|F| > 1, p < 0.01")) +
  labs(y = "log(fold change)", x = "Mean normalized counts") + 
  theme_bw() + 
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "top")
ggsave(filename = paste0("res/", format(Sys.Date(), "%Y-%m-%d"), "_MAplot.tiff"), device = "tiff", width = 11, height = 8, dpi = 600)

res_ic18 <- results(dds,
                   contrast = contrast_ic18,
                   pAdjustMethod = "fdr",
                   cooksCutoff = FALSE)

res_ic18 <- res_ic18 %>%
  as_tibble() %>%
  tibble::add_column(., gene = gkey$SYMBOL, .before = 1) %>% 
  tibble::add_column(., entrez_id = gkey$ENTREZID, .before = 2)

# viral specific genes
fold_thr <- 1
pval <- 0.05

# virus_sig <- with(res_v18, which(padj < 0.05 & log2FoldChange > 1))
virus_sig <- res_v18 %>% 
  dplyr::filter(., abs(log2FoldChange) > fold_thr, padj < pval) 

ic_sig <- res_ic18 %>% 
  dplyr::filter(., abs(log2FoldChange) > fold_thr, padj < pval) 


# virus_specific <- setdiff(virus_sig, ic_sig)
virus_specific <- virus_sig %>% 
  tibble::add_column(., ic_gene = 0) %>% 
  dplyr::mutate(., ic_gene = replace(ic_gene, which(gene %in% ic_sig$gene), 1)) %>% 
  dplyr::filter(., ic_gene == 0)

res_v18 %>% 
  dplyr::mutate(sig, sig = replace(sig, which(entrez_id %in% virus_specific$entrez_id), 2)) %>%
  ggplot() + 
  geom_point(aes(y = log2FoldChange, x = log10(baseMean), color = as.factor(sig)), alpha = 0.5, size = 3) + 
  ylim(c(-max(abs(res_v18$log2FoldChange)), max(abs(res_v18$log2FoldChange)))) +
  scale_color_viridis_d(breaks = c(0, 1, 2), name = NULL, labels = c("Not significant", "Significant, general", "Significant, virus specific")) +
  # scale_color_manual(breaks = c(0, 1, 2), name = NULL, values = c("black", "red", "tan"), labels = c("Not significant", "|F| > 1, p < 0.01", "Significant, virus specific")) +
  labs(y = "log(fold change)", x = "Mean normalized counts") + 
  theme_bw() + 
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "top")
ggsave(filename = filename = paste0("res/", Sys.Date(),"_virus-specific_MAplot.", dev),
       device = dev, 
       width = 11, 
       height = 8, 
       dpi = 600)

# build gene set for DRUID ----
druid_genes <- cbind(virus_specific$log2FoldChange,
                     virus_specific$padj)

flu_druid_2 <- DRUID::concoct(dge_matrix = druid_genes, 
                            num_random = 10000, 
                            druid_direction = "neg", 
                            fold_thr = 0, 
                            pvalue_thr = 0.05, 
                            entrez = virus_specific$entrez_id)

dev <- "tiff"
flu_druid %>% 
  dplyr::filter(., druid_score > 5, number_matches > 1) %>% 
  # dplyr::slice(1:1000) %>%
  ggplot() + 
  geom_point(aes(x = druid_score, y = drug_name, color = druid_score), alpha = 0.5, size = 4) +
  scale_color_viridis_c() +
  # facet_grid(. ~ data_source, scales = "free") +
  labs(x = "DRUID score", y = "Drug") + 
  theme_bw() + 
  theme(axis.title = element_text(size = 24, color = "black"),
        strip.text.x = element_text(size = 12, color = "black"),
        axis.text = element_text(size = 12, color = "black"),
        axis.ticks.length = unit(0.2, "cm"))
ggsave(filename = paste0("res/", Sys.Date(),"_DRUID-chip-data.", dev),
       device = dev,
       height = 8,
       width = 11,
       dpi = 600)

flu_grove <- grove::stir(dge_matrix = druid_genes, 
                         effect = "neg", 
                         fold_thr = 1, 
                         pvalue_thr = 0.05, 
                         entrez = gkey$ENTREZID, 
                         num_random = 10000, 
                         num_combos = 10)

# tissue specificity using stouffer's method
#
source("R/lung_specificity.R")

# pathway enrichment ----
# using oranges (fisher exact test implementation)
library(oranges)

chip_enr_ora <- oranges::oranges(query_entrez = virus_specific$entrez_id, universe_entrez = gkey$ENTREZID)
dev <- "tiff"
chip_enr_ora %>% 
  dplyr::filter(., data_source == "go:bp") %>% 
  dplyr::arrange(., desc(-log10(padj))) %>% 
  dplyr::slice(1:25) %>%
  ggplot() + 
  geom_point(aes(x = -log10(padj), y = forcats::fct_reorder(name, -log10(padj)), color = -log10(padj), size = number_genes), alpha = 0.5) +
  scale_color_viridis_c() + 
  labs(x = "-log10(adjusted p-value)", y = "Pathway name") +
  theme_bw() + 
  theme(axis.title = element_text(size = 24, color = "black"),
        axis.text = element_text(size = 12, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none")
ggsave(filename = paste0("res/", Sys.Date(),"_chip-data-ORANGES-enrichment.", dev),
       device = dev,
       height = 8,
       width = 11,
       dpi = 600)

# 
flu_enr_rpegeos <- rpegeos::enrich_geneset(gene_set = cbind(virus_specific$entrez_id, virus_specific$log2FoldChange, virus_specific$padj))
dev <- "tiff"
chip_enr %>%
  dplyr::slice(1:20) %>%
  ggplot() +
  geom_point(aes(y = enrichment_score, x = forcats::fct_reorder(geneset, enrichment_score), size = number_genes), alpha = 0.5, color = "#ff6600") +
  labs(x = "Gene set", y = "Enrichment score") +
  coord_flip() + 
  theme_bw() + 
  theme(axis.title = element_text(size = 12, color = "black"),
        axis.text = element_text(size = 12, color = "black"),
        legend.position = "none")
ggsave(filename = paste0("res/", Sys.Date(),"_chip-data-enrichment.", dev),
       device = dev,
       height = 8,
       width = 11,
       dpi = 600)


flu_enr_cp <- clusterProfiler::enrichGO(gene = virus_specific$entrez_id, ont = "BP", pAdjustMethod = "fdr", universe = gkey$ENTREZID, OrgDb = org.Hs.eg.db)

# expression of specific genes ----
# neu3 and b2m
tmp_df <- tibble(gene_name = rep(gkey$SYMBOL, ncol(count_data)),
                 gene_entrez = rep(gkey$ENTREZID, ncol(count_data)),
                 counts = as.vector(as.matrix(count_data)),
                 group = as.vector(sapply(sample_data$group, rep, length(gene_ids))),
                 timepoint = as.vector(sapply(sample_data$time_point, rep, length(gene_ids))),
                 group_time = as.vector(sapply(sample_data$group_time, rep, length(gene_ids))))

dev <- "tiff"
tmp_df %>% 
  dplyr::filter(., timepoint == 18) %>%
  dplyr::filter(., gene_name %in% intersect(res_v18$gene[which(res_v18$sig == 1)], lung_genes)) %>% 
  ggplot() + 
  geom_boxplot(aes(x = group, y = counts, fill = group)) + 
  geom_jitter(aes(x = group, y = counts), width = 0.1) + 
  scale_fill_viridis_d() + 
  facet_grid(. ~ gene_name, scales = "free") +
  theme_bw() + 
  theme(axis.text = element_text(size = 18, color = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 36, color = "black"), 
        axis.ticks.length = unit(0.2, "cm"),
        strip.text.x = element_text(size = 12, color = "black"),
        panel.grid = element_blank())
ggsave(filename = paste0("res/", Sys.Date(),"_LUNG-SPECIFIC_expression.", dev),
       device = dev,
       height = 8,
       width = 11,
       dpi = 600)



dev <- "tiff"
tmp_df %>% 
  dplyr::filter(., timepoint == 18) %>%
  dplyr::filter(., gene_name %in% virus_biomarkers$gene[virus_biomarkers$secreted == 1]) %>% 
  ggplot() + 
  geom_boxplot(aes(x = group, y = counts, fill = group)) + 
  geom_jitter(aes(x = group, y = counts), width = 0.1) + 
  scale_fill_viridis_d() + 
  facet_wrap(. ~ gene_name, scales = "free_y") +
  theme_bw() + 
  theme(axis.text = element_text(size = 18, color = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 36, color = "black"), 
        axis.ticks.length = unit(0.2, "cm"),
        strip.text.x = element_text(size = 12, color = "black"),
        panel.grid = element_blank())
ggsave(filename = paste0("res/", Sys.Date(),"_SECRETED-BIOMARKERS_expression.", dev),
       device = dev,
       height = 8,
       width = 11,
       dpi = 600)

# genes of interest
gint <- unique(c(virus_biomarkers %>% 
            dplyr::filter(., secreted == 1 | signal_peptide == 1, padj < 0.05, abs(log2FoldChange) > 1) %>% 
            dplyr::select(., gene) %>% 
            as.matrix %>% 
            as.vector,
          intersect(lung_genes, res_v18$gene[res_v18$padj < 0.05 & abs(res_v18$log2FoldChange) > 1])))
dev <- "tiff"
tmp_df %>% 
  dplyr::filter(., timepoint == 18) %>%
  dplyr::filter(., gene_name %in% gint) %>% 
  ggplot() + 
  geom_boxplot(aes(x = group, y = counts, fill = group)) + 
  geom_jitter(aes(x = group, y = counts), width = 0.1) + 
  scale_fill_viridis_d() + 
  facet_wrap(. ~ gene_name, scales = "free_y") +
  theme_bw() + 
  theme(axis.text = element_text(size = 18, color = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 36, color = "black"), 
        axis.ticks.length = unit(0.2, "cm"),
        strip.text.x = element_text(size = 12, color = "black"),
        panel.grid = element_blank())
ggsave(filename = paste0("res/", Sys.Date(),"_GENES-INTEREST_expression.", dev),
       device = dev,
       height = 8,
       width = 11,
       dpi = 600)