# pagerank idea

load("data/2018-04-30_ttrust_data.RData")

net <- ttrust_data %>% 
  dplyr::mutate(., tf2 = tf) %>% 
  dplyr::mutate(., count_refs = sapply(reference, 
                                       function(y) length(strsplit(y, ";")[[1]]))) %>%
  dplyr::select(., target, tf2, count_refs) %>% 
  igraph::graph_from_data_frame(., directed = TRUE)

# pagerank of genes ----
tmp1 <- igraph::page_rank(graph = net, directed = TRUE)

# genes we're interested in ----
xx <- sort(tmp1$vector[match(virus_specific$gene, 
                             names(tmp1$vector), 
                             nomatch = FALSE)], 
           decreasing = TRUE)

net2 <- igraph::as_data_frame(net)
net2 %>% 
  dplyr::filter(., from %in% names(xx) | to %in% names(xx))


# annotate genes ----
yy <- AnnotationDbi::select(x = org.Hs.eg.db, keys = names(xx), columns = "ENTREZID", keytype = "SYMBOL")


# get expression data for those genes ----
z1 <- match(yy$SYMBOL, gene_ids)
z1 <- z1[!is.na(z1)]
z2 <- cbind(gene_ids[z1], res_v18$log2FoldChange[z1], res_v18$padj[z1])
z3 <- AnnotationDbi::select(x = org.Hs.eg.db, 
                            keys = z2[, 1],
                            keytype = "SYMBOL",
                            columns = "ENTREZID")
z4 <- cbind(z3$ENTREZID, z2[, c(2, 3)])
z4 <- apply(z4, 2, as.numeric)

# enrich ----
# z5 <- rpegeos::enrich_geneset(z4)
z5 <- oranges::oranges(query_entrez = z4[, 1], universe_entrez = as.vector(AnnotationDbi::mapIds(x = org.Hs.eg.db, keys = names(tmp1$vector), column = "ENTREZID", keytype = "SYMBOL")))

dev <- "tiff"
z5 %>% 
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

a1 <- ttrust_data %>%
  dplyr::filter(., target %in% virus_specific$gene | tf %in% virus_specific$gene)

net <- a1 %>% 
  igraph::graph_from_data_frame()

igraph::vertex_attr(net) <- list(name = unique(c(a1$tf, a1$target)),
                                 color = rep("red", length(igraph::V(net))), 
                                 size = rep(4, length(igraph::V(net))))

plot(net)

