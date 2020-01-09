a1 <- flu_enr_ora
a3 <- chip_enr_ora

a1 <- a1 %>% 
  tibble::add_column(., group = "zaas data") %>% 
  tibble::add_column(., rank = seq(1, nrow(a1))) %>% 
  dplyr::filter(., data_source == "go:bp") %>%
  dplyr::arrange(., desc(-log10(padj))) %>% 
  dplyr::slice(1:25)

a3 <- a3 %>% 
  tibble::add_column(., group = "transwell data") %>% 
  tibble::add_column(., rank = seq(1, nrow(a3))) %>%
  dplyr::filter(., data_source == "go:bp") %>%
  dplyr::arrange(., desc(-log10(padj))) %>% 
  dplyr::slice(1:25)


a4 <- rbind(a1, a3)

dev <- "tiff"
a4 %>% 
  # dplyr::filter(., group != "network data") %>%
  # dplyr::filter(., data_source == "go:bp") %>%
  ggplot() + 
  geom_point(aes(y = forcats::fct_reorder(name, -log10(padj)), x = -log10(padj), color = group, size = number_genes)) +
  scale_color_manual(values = c("#ff9933", "#cc99cc", "#99cc99")) +
  labs(x = "Gene set", y = "Enrichment score") +
  # coord_flip() + 
  facet_grid(. ~ group, scales = "free") +
  theme_bw() + 
  theme(axis.title = element_text(size = 24, color = "black"),
        axis.text = element_text(size = 12, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none")
ggsave(filename = paste0("res/", Sys.Date(),"_combine-enrichments.", dev),
       device = dev,
       height = 8,
       width = 11,
       dpi = 600)


# plot virus specific + lung specific genes ----
b1 <- intersect(virus_specific$gene, lung_genes)
b2 <- res_v18$log2FoldChange[which(gene_ids %in% b1)] # <-- chip data
b3 <- y1$logFC[which(gene_mappings_copd$SYMBOL %in% b1)] # <-- copd data
b4 <- xx1$logFC[which(G$SYMBOL %in% b1)] # <-- zaas data, flu
b5 <- xx2$logFC[which(G$SYMBOL %in% b1)] # <-- zaas data, rsv
b6 <- xx3$logFC[which(G$SYMBOL %in% b1)] # <-- zaas data, hrv

dge_groups <- tibble(group = c(rep("this study", length(b2)),
                               rep("copd data (benam et al)", length(b3)),
                               rep("flu infection data (zaas et al)", length(b4)),
                               rep("rsv infection data (zaas et al)", length(b5)),
                               rep("hrv infection data (zaas et al)", length(b6))),
                     gene_name = c(gene_ids[which(gene_ids %in% b1)], 
                                   gene_mappings_copd$SYMBOL[which(gene_mappings_copd$SYMBOL %in% b1)],
                                   G$SYMBOL[which(G$SYMBOL %in% b1)],
                                   G$SYMBOL[which(G$SYMBOL %in% b1)],
                                   G$SYMBOL[which(G$SYMBOL %in% b1)]),
                     fold_data = c(b2, b3, b4, b5, b6))

dge_groups %>% 
  dplyr::filter(., group != "copd data (benam et al)", group != "rsv infection data (zaas et al)", group != "hrv infection data (zaas et al)") %>%
  ggplot() + 
  geom_point(aes(x = group, y = fold_data, color = group), size = 5, alpha = 0.5) + 
  scale_color_manual(values = c("purple", "orange")) +
  facet_wrap(. ~ gene_name) + 
  labs(x = NULL, y = "log2(fold difference to control)") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 12, color = "black", angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        legend.position = "none")
ggsave(filename = paste0("res/", Sys.Date(),"_LUNG-SPECIFIC-COMPARISON_fold_changes.png"),
       device = "png",
       height = 8,
       width = 11,
       dpi = 600)


# secreted virus specific ----
d2 <- res_v18$log2FoldChange[which(gene_ids %in% virus_biom_sec)] # <-- chip data
d4 <- xx1$logFC[which(G$SYMBOL %in% virus_biom_sec)] # <-- zaas data, flu

sec_genes <- tibble(group = c(rep("this study", length(d2)),
                               rep("flu infection data (zaas et al)", length(d4))),
                     gene_name = c(gene_ids[which(gene_ids %in% virus_biom_sec)], 
                                   G$SYMBOL[which(G$SYMBOL %in% virus_biom_sec)]),
                     fold_data = c(d2, d4))

sec_genes %>% 
  ggplot() + 
  geom_point(aes(x = forcats::fct_reorder(gene_name, fold_data, .desc = TRUE), y = fold_data, color = group), size = 3, alpha = 0.5) + 
  geom_hline(yintercept = 1, color = "black", lty = 2) + 
  scale_color_manual(values = c("purple", "orange")) +
  labs(x = "Secreted biomarker", y = "log2(fold difference to control)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12))
ggsave(filename = paste0("res/", Sys.Date(),"_SECRETED-BIOMARKERS_comparison.png"),
       device = "png",
       height = 8,
       width = 11,
       dpi = 600)



