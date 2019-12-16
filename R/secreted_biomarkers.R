# signal peptides ----
sig_pep <- read_xlsx("~/work/data/2019-03-07_signal_peptides_human.xlsx")

sig_peps_names <- AnnotationDbi::select(x = org.Hs.eg.db, keys = as.character(sig_pep$`Accession Number`), keytype = "UNIPROT", columns = c("ENTREZID", "SYMBOL"))
sig_peps_names <- sig_peps_names[which(!is.na(sig_peps_names$ENTREZID)), ]

# differentially abundant signal peptides in our data
# virus_sig_peps <- intersect(sig_peps_names$SYMBOL, gene_ids[virus_specific])


# secreted proteins ----
# based on HPA predicted data
hpa_pred <- read_delim("~/work/data/2019-03-07_HPA_predicted_secreted_proteins.tsv", delim = "\t")
hpa_sec <- hpa_pred[grep("secreted", hpa_pred$`Protein class`), ]

# virus_secp <- intersect(hpa_pred$Gene, gene_ids[virus_specific])

# predicted secreted biomarkers ----
# virus_biom_sec <- union(virus_secp, virus_sig_peps)
virus_biomarkers <- virus_specific %>% 
  tibble::add_column(., signal_peptide = 0) %>% 
  tibble::add_column(., secreted = 0) %>%
  dplyr::mutate(., signal_peptide = replace(signal_peptide, which(gene %in% sig_peps_names$SYMBOL), 1)) %>% 
  dplyr::mutate(., secreted = replace(secreted, which(gene %in% hpa_pred$Gene), 1))

# plot ----
dev <- "tiff"
virus_biomarkers %>% 
  dplyr::filter(., signal_peptide == 1 | secreted == 1) %>%
  ggplot() + 
  geom_point(aes(x = forcats::fct_reorder(gene, log2FoldChange, .desc = TRUE), y = log2FoldChange, color = log2FoldChange, size = -log10(padj)), alpha = 0.5) + 
  scale_color_viridis_c() + 
  labs(x = NULL, y = "log2(Fold change)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 24, color = "black"))
ggsave(filename = paste0("res/", Sys.Date(),"_virus-specific-biomarkers.", dev),
       device = dev, 
       width = 11, 
       height = 8, 
       dpi = 600)



