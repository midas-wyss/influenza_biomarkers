# get data
load("res/2019-04-09_diffnet_zaas.RData")
load(file = "data/2018-04-30_ttrust_data.RData")

# filtering
tmp2 <- DN %>% 
  dplyr::mutate(., x_name = AnnotationDbi::select(x = org.Hs.eg.db, keys = as.character(x), keytype = "ENTREZID", columns = "SYMBOL")$SYMBOL) %>%
  dplyr::mutate(., y_name = AnnotationDbi::select(x = org.Hs.eg.db, keys = as.character(y), keytype = "ENTREZID", columns = "SYMBOL")$SYMBOL) %>%
  dplyr::filter(., !is.na(x_name), !is.na(y_name)) %>% 
  dplyr::mutate(., x_tf = x_name %in% ttrust_data$tf) %>%
  dplyr::mutate(., x_tf = replace(x_tf, x_tf == TRUE, 1), x_tf = replace(x_tf, x_tf == FALSE, 0)) %>%
  dplyr::mutate(., y_tf = y_name %in% ttrust_data$tf) %>%
  dplyr::mutate(., y_tf = replace(y_tf, y_tf == TRUE, 1), y_tf = replace(y_tf, y_tf == FALSE, 0))


# look at where the virus specific genes and lung specific genes are
xx <- gene_ids[virus_specific]
xx <- AnnotationDbi::select(x = org.Hs.eg.db, keys = xx, columns = "ENTREZID", keytype = "SYMBOL")

tmp <- tmp2 %>% 
  dplyr::filter(., x %in% gkey$ENTREZID | y %in% gkey$ENTREZID) %>% 
  dplyr::filter(., y %in% gkey$ENTREZID | x %in% gkey$ENTREZID) %>% 
  dplyr::filter(., change_type != "edge present, same sign") %>% 
  dplyr::mutate(., x_lung = x_name %in% lung_genes) %>%
  dplyr::mutate(., x_lung = replace(x_lung, x_lung == TRUE, 1), x_lung = replace(x_lung, x_lung == FALSE, 0)) %>%
  dplyr::mutate(., y_lung = y_name %in% lung_genes) %>%
  dplyr::mutate(., y_lung = replace(y_lung, y_lung == TRUE, 1), y_lung = replace(y_lung, y_lung == FALSE, 0)) %>% 
  # dplyr::filter(., x_lung == 1 | y_lung == 1) %>% 
  dplyr::filter(., x_name %in% gkey$SYMBOL, y_name %in% gkey$SYMBOL) %>%
  dplyr::mutate(., x_virus = x_name %in% virus_specific$gene) %>%
  dplyr::mutate(., x_virus = replace(x_virus, x_virus == TRUE, 1), x_virus = replace(x_virus, x_virus == FALSE, 0)) %>%
  dplyr::mutate(., y_virus = y_name %in% virus_specific$gene) %>%
  dplyr::mutate(., y_virus = replace(y_virus, y_virus == TRUE, 1), y_virus = replace(y_virus, y_virus == FALSE, 0))


tmp3 <- tmp %>%
  dplyr::mutate(., viral_gene = x_virus + y_virus) %>%
  dplyr::mutate(., tf_gene = x_tf + y_tf) %>%
  dplyr::mutate(., lung_gene = x_lung + y_lung) %>%
  dplyr::filter(., viral_gene != 0, tf_gene != 0) %>% 
  dplyr::select(., -p_val, -change_type, - edge_score, -x, -y, -x_virus, -y_virus)

sw <- which(tmp3$y_tf == 1 & tmp3$x_tf == 0)
tmp4 <- tmp3
tmp4$x_name[sw] <- tmp3$y_name[sw]
tmp4$y_name[sw] <- tmp3$x_name[sw]

tmp4 %>%  
  dplyr::group_by(., lung_gene, x_name) %>% 
  dplyr::count() %>% 
  dplyr::filter(., n > 2) %>%
  ggplot() +
  geom_point(aes(x = x_name, y = n, color = as.factor(lung_gene)), size = 4, alpha = 0.5) + 
  scale_color_viridis_d(breaks = c(0, 1), labels = c("Not lung specific", "Lung specific"), name = "Tissue Specificity") +
  labs(x = NULL, y = NULL) +
  # facet_grid(. ~ lung_gene, scales = "free") +
  theme_bw() +
  theme(axis.ticks.length = unit(0.2, "cm"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"))

  
tmp4 %>%  
  dplyr::group_by(., lung_gene, x_name) %>% 
  dplyr::count() %>% 
  dplyr::filter(., n > 2) %>%
  dplyr::mutate(., fold = res_v18$log2FoldChange[which(res_v18$gene %in% x_name)]) %>%
  ggplot() +
  geom_point(aes(x = forcats::fct_reorder(x_name, fold, .desc = TRUE), y = n, color = as.factor(lung_gene)), size = 4, alpha = 0.5) + 
  scale_color_viridis_d(breaks = c(0, 1), labels = c("Not lung specific", "Lung specific"), name = "Tissue Specificity") +
  labs(x = NULL, y = NULL) +
  # facet_grid(. ~ lung_gene, scales = "free") +
  theme_bw() +
  theme(axis.ticks.length = unit(0.2, "cm"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"))


node_names <- c(tmp$x, tmp$y)
node_names <- unique(node_names)
node_names <- AnnotationDbi::select(x = org.Hs.eg.db, keys = as.character(node_names), columns = "SYMBOL", keytype = "ENTREZID")

node_names <- node_names %>% 
  dplyr::mutate(., lung_gene = node_names$SYMBOL %in% lung_genes) %>%
  dplyr::mutate(., lung_gene = replace(lung_gene, lung_gene == TRUE, 1), lung_gene = replace(lung_gene, lung_gene == FALSE, 0)) %>%
  dplyr::mutate(., virus_gene = node_names$SYMBOL %in% gene_ids[virus_specific]) %>%
  dplyr::mutate(., virus_gene = replace(virus_gene, virus_gene == TRUE, 1), virus_gene = replace(virus_gene, virus_gene == FALSE, 0)) %>%
  dplyr::mutate(., log2fc = res_v18$log2FoldChange[match(node_names$SYMBOL, gene_ids)]) %>%
  dplyr::mutate(., p_val = res_v18$padj[match(node_names$SYMBOL, gene_ids)])



# write files to load into cytoscape ----
# maybe do something for igraph instead? 
write_delim(x = tmp, "res/2019-04-08_diffnet_results.txt", delim = "\t", col_names = TRUE)
write_delim(x = node_names, "res/2019-04-08_diffnet_results_nodenames.txt", delim = "\t", col_names = TRUE)


# enrichment of nodes in subnetwork ----
subnet_gset <- cbind(node_names$ENTREZID, node_names$log2fc, node_names$p_val)
subnet_gset <- apply(subnet_gset, 2, as.numeric)
subnet_enr <- rpegeos::enrich_geneset(gene_set = subnet_gset)