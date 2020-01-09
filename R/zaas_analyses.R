# packages ----
library(hgu133plus2.db)
library(limma)
library(tidyverse)

# functions ----
source("/Volumes/HOME/scripts/r/expression_analysis/limma_dge.R")


# DREAM respiratory challenge data ----
load("data/public_data/zaas_data_05152018.RData")
load("data/public_data/2019-02-12_GSE30273.RData")

E <- data[[1]]
G <- data[[2]]
S <- data[[3]]

# separate viral groups ----
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

# differential expression ----
xx1 <- limma_dge(expression_data = E,
                 ctrIds = flu_peak_asymp,
                 caseIds = flu_peak_symp)

xx2 <- limma_dge(expression_data = E,
                 ctrIds = rsv_peak_asymp,
                 caseIds = rsv_peak_symp)

xx3 <- limma_dge(expression_data = E,
                 ctrIds = hrv_peak_asymp,
                 caseIds = hrv_peak_symp)


diff_genes_flu <- which(xx1$adj.P.Val < 0.01 & abs(xx1$logFC) > 1)
diff_genes_rsv <- which(xx2$adj.P.Val < 0.01 & abs(xx2$logFC) > 1)
diff_genes_hrv <- which(xx3$adj.P.Val < 0.01 & abs(xx3$logFC) > 1)


# flu biomarkers ----
# these are genes that are differentially expressed in flu but not
# in the other 2 infection models
diff_genes_flu <- setdiff(diff_genes_flu,
                          union(diff_genes_hrv, diff_genes_rsv))

diff_res <- tibble::tibble(gene_id=G$ENTREZID[diff_genes_flu],
                       gene_symbol=G$SYMBOL[diff_genes_flu],
                       fold_change=xx1$logFC[diff_genes_flu],
                       fdr_pvalue=xx1$adj.P.Val[diff_genes_flu])

# gene set enrichment with rpegeos ----
library(rpegeos)
library(oranges)
gset <- cbind(G$ENTREZID[diff_genes_flu], xx1$logFC[diff_genes_flu], xx1$adj.P.Val[diff_genes_flu])
gset <- apply(gset, 2, as.numeric)
flu_enr <- rpegeos::enrich_geneset(gene_set = gset)

flu_enr_ora <- oranges::oranges(query_entrez = G$ENTREZID[diff_genes_flu], universe_entrez = G$ENTREZID)

dev <- "tiff"
flu_enr_ora %>% 
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
ggsave(filename = paste0("res/", Sys.Date(),"_zaas-data-ORANGES-enrichment.", dev),
       device = dev,
       height = 8,
       width = 11,
       dpi = 600)


# biomarker discovery ----
# compare biomarkers here with ones identified on chip data from study


# DRUID ----
zaas_genes <- cbind(xx1$logFC, xx1$adj.P.Val)

zaas_druid <- DRUID::concoct(dge_matrix = zaas_genes, 
                            num_random = 10000, 
                            druid_direction = "neg", 
                            fold_thr = 1, 
                            pvalue_thr = 0.01, 
                            entrez = G$ENTREZID)

dev <- "tiff"
zaas_druid %>% 
  dplyr::filter(., druid_score > 5.5, number_matches > 1) %>% 
  # dplyr::slice(1:1000) %>%
  ggplot() + 
  geom_point(aes(x = druid_score, y = drug_name, color = druid_score, size = number_matches), alpha = 0.5) +
  scale_color_viridis_c() +
  facet_grid(. ~ data_source, scales = "free") +
  labs(x = "DRUID score", y = "Drug") + 
  theme_bw() + 
  theme(axis.title = element_text(size = 24, color = "black"),
        strip.text.x = element_text(size = 12, color = "black"),
        axis.text = element_text(size = 12, color = "black"),
        axis.ticks.length = unit(0.2, "cm"))
ggsave(filename = paste0("res/", Sys.Date(),"_DRUID-zaas-data.", dev),
       device = dev,
       height = 8,
       width = 11,
       dpi = 600)

x1 <- zaas_druid %>% dplyr::filter(., druid_score > 5) %>% dplyr::select(drug_name) %>% unique %>% as.matrix %>% as.vector
x2 <- flu_druid %>% dplyr::filter(., druid_score > 5) %>% dplyr::select(drug_name) %>% unique %>% as.matrix %>% as.vector