# run DRUID on chip signature
# with Alberto de la Fuente and Ratnakar Potla

# libraries
library(readxl)
library(org.Hs.eg.db)
library(DRUID)


# load data
deg_list <- read_xlsx(path = "data/NEW_UNIQUE_LIST.xlsx", sheet = 1)

# map genes
genes <- AnnotationDbi::select(x = org.Hs.eg.db, keys = as.character(deg_list$Entrez), keytype = "ENTREZID", columns = c("SYMBOL", "GENENAME"))

# generate gene set
gs <- cbind(deg_list$Log2FC, rep(0, nrow(deg_list)))

# run DRUID
# alf_druid <- run_conddr(conddr_data = "/Volumes/HOME/scripts/r/conddr/data/conddr_data_ndrugs=8422_neffects=41412_saveDate=03082018.RData",
#                            dge_data = gs,
#                            pvalue_thr = 0.01,
#                            log2_fold_thr = 0,
#                            gene_symbols = genes$SYMBOL,
#                            gene_entrez = genes$ENTREZID,
#                            num_random = 10000,
#                            match_effect = "neg")
# 
flu_druid <- concoct(dge_matrix = gs, 
                     tfidf_matrix = cmap_druid$tfidf, 
                     tfidf_crossproduct = cmap_druid$cpm, 
                     num_random = 10000, 
                     druid_direction = "neg", 
                     fold_thr = 0, 
                     pvalue_thr = 0.01, 
                     entrez = genes$ENTREZID)

flu_druid <- flu_druid %>% 
  tibble::add_column(., drug_name = cmap_druid$drugs$name, .before = 1) %>%
  tibble::add_column(., concentration = cmap_druid$drugs$concentration, .before = 2) %>%
  tibble::add_column(., cell_line = cmap_druid$drugs$cell_line, .before = 3)



# genes --> celular compartment
g2 <- AnnotationDbi::select(x = org.Hs.eg.db, keys = as.character(deg_list$Entrez), keytype = "ENTREZID", columns = c("SYMBOL", "GENENAME", "GO"))

x1 <- which(g2$ONTOLOGY == "CC")
x2 <- unique(g2$GO[x1])

g2cc <- vector(mode = "list", length = length(x2))

for (i in seq(1, length(x2))) {
  g2cc[[i]] <- data_frame(go_term = x2[i],
                          number_genes = length(unique(g2$ENTREZID[which(g2$GO == x2[i])])))
}
g2cc <- dplyr::bind_rows(g2cc)

g2cc %>% 
  dplyr::filter(., number_genes > 10) %>% 
  ggplot() +
  geom_bar(aes(x = fct_reorder(go_term, desc(number_genes)), y = number_genes), position = "dodge", color = "black", stat = "identity") + 
  theme_bw()
