# pathway revealer ----
# given a gene set, find pathways that those genes are involved in and color the pathways appropriately based on gene expression data

# kegg pathways (from rpegeos)
load(file = "~/work/data/pathway_sets.RData")

#
#
gint <- c("GHRL", "MMP1", "MMP8")
x <- AnnotationDbi::select(x = org.Hs.eg.db, keys = gint, keytype = "SYMBOL", columns = "ENTREZID")

#
# 
# map genes to pathway matrix columns
pmat <- pathway_sets$membership_matrix
pmat <- pmat[which(pathway_sets$pathway_info$pathway_source == "kegg"), ]
pnames <- pathway_sets$pathway_info[which(pathway_sets$pathway_info$pathway_source == "kegg"), ]
gids <- which(colnames(pmat) %in% x$ENTREZID)

pints <- vector(mode = "list", length = length(gids))
for (i in seq(1, length(gids))) {
  a1 <- pmat[, gids[i]]
  a2 <- which(a1 == 1)
  if (length(a2) == 1) {
    pints[[i]] <- a2
  } else {
    a3 <- pmat[a2, ]
    a4 <- Matrix::rowSums(a3)
    a5 <- which(a4 < 100)
    pints[[i]] <- a2[a5]
  }
}
pints <- unique(unlist(pints))

sub_pmat <- pmat[pints, ]
sub_names <- pnames[pints, ]
sub_names <- as.vector(as.matrix(sub_names[, 2]))
nix <- c(1, 3)
sub_pmat <- sub_pmat[-nix, ]
sub_names <- sub_names[-nix]

nix <- which(Matrix::colSums(sub_pmat) == 0)
if(length(nix) != 0) {
  sub_pmat <- sub_pmat[, -nix]  
}

# paint matrix
paint_mat <- vector(mode = "list", length = nrow(sub_pmat))
for (i in seq(1, nrow(sub_pmat))) {
  y1 <- names(which(sub_pmat[i, ] == 1))
  z1 <- AnnotationDbi::select(x = org.Hs.eg.db, 
                              keys = y1, 
                              keytype = "ENTREZID", 
                              columns = "SYMBOL")
  y2 <- match(y1, gkey$ENTREZID)
  y3 <- cbind(z1, res_v18$log2FoldChange[y2], res_v18$padj[y2])
  colnames(y3) <- c("entrez_id", "symbol", "logFC", "padj")
  y3 <- as_tibble(y3)
  y3 <- y3 %>%  
    tibble::add_column(., pathway = sub_names[i], .before = 1) %>%
    tibble::add_column(., direction = 0) %>% 
    dplyr::mutate(., direction = replace(direction, logFC > 0 & padj < 0.05, 1)) %>%
    dplyr::mutate(., direction = replace(direction, logFC < 0 & padj < 0.05, -1))
  paint_mat[[i]] <- y3
}
paint_mat <- dplyr::bind_rows(paint_mat)


paint_mat %>% 
  # dplyr::filter(., pathway == "ppar_signaling_pathway") %>%
  dplyr::filter(., !is.na(padj)) %>%
  ggplot() + 
  geom_tile(aes(x = symbol, y = pathway, fill = logFC), color = "black") +
  # geom_point(aes(x = forcats::fct_reorder(symbol, logFC), y = pathway, fill = logFC), fill = "black") + 
  # scale_color_gradient2(low = "darkblue", mid = "white", high = "red", midpoint = 0, na.value = "white") + 
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "red", midpoint = 0, na.value = "white") + 
  coord_fixed(ratio = 6) +
  labs(x = "Gene", y = NULL) +
  # facet_grid(. ~ pathway, scales = "free") +
  theme_bw() + 
  theme(axis.text = element_text(size = 12, color = "black"),
        # axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.x = element_blank(),
        panel.grid = element_blank())

##################
# MY BIOMARKERS
#
#
gint <- c("WARS", "TNFSF10")
x <- AnnotationDbi::select(x = org.Hs.eg.db, keys = gint, keytype = "SYMBOL", columns = "ENTREZID")

#
# 
# map genes to pathway matrix columns
pmat <- pathway_sets$membership_matrix
pmat <- pmat[which(pathway_sets$pathway_info$pathway_source == "kegg"), ]
pnames <- pathway_sets$pathway_info[which(pathway_sets$pathway_info$pathway_source == "kegg"), ]
gids <- which(colnames(pmat) %in% x$ENTREZID)

pints <- vector(mode = "list", length = length(gids))
for (i in seq(1, length(gids))) {
  a1 <- pmat[, gids[i]]
  a2 <- which(a1 == 1)
  if (length(a2) == 1) {
    pints[[i]] <- a2
  } else {
    a3 <- pmat[a2, ]
    a4 <- Matrix::rowSums(a3)
    a5 <- which(a4 < 100)
    pints[[i]] <- a2[a5]
  }
}
pints <- unique(unlist(pints))

sub_pmat <- pmat[pints, ]
sub_names <- pnames[pints, ]
sub_names <- as.vector(as.matrix(sub_names[, 2]))
nix <- 3
sub_pmat <- sub_pmat[-nix, ]
sub_names <- sub_names[-nix]

nix <- which(Matrix::colSums(sub_pmat) == 0)
if(length(nix) != 0) {
  sub_pmat <- sub_pmat[, -nix]  
}

# paint matrix
paint_mat <- vector(mode = "list", length = nrow(sub_pmat))
for (i in seq(1, nrow(sub_pmat))) {
  y1 <- names(which(sub_pmat[i, ] == 1))
  z1 <- AnnotationDbi::select(x = org.Hs.eg.db, 
                              keys = y1, 
                              keytype = "ENTREZID", 
                              columns = "SYMBOL")
  y2 <- match(y1, gkey$ENTREZID)
  y3 <- cbind(z1, res_v18$log2FoldChange[y2], res_v18$padj[y2])
  colnames(y3) <- c("entrez_id", "symbol", "logFC", "padj")
  y3 <- as_tibble(y3)
  y3 <- y3 %>%  
    tibble::add_column(., pathway = sub_names[i], .before = 1) %>%
    tibble::add_column(., direction = 0) %>% 
    dplyr::mutate(., direction = replace(direction, logFC > 0 & padj < 0.05, 1)) %>%
    dplyr::mutate(., direction = replace(direction, logFC < 0 & padj < 0.05, -1))
  paint_mat[[i]] <- y3
}
paint_mat <- dplyr::bind_rows(paint_mat)


paint_mat %>% 
  # dplyr::filter(., pathway == "ppar_signaling_pathway") %>%
  dplyr::filter(., !is.na(padj)) %>%
  ggplot() + 
  # geom_tile(aes(x = symbol, y = pathway, fill = logFC), color = "black") +
  geom_tile(aes(x = forcats::fct_reorder(symbol, logFC), y = pathway, fill = logFC), color = "black") +
  # geom_point(aes(x = forcats::fct_reorder(symbol, logFC), y = pathway, fill = logFC), fill = "black") + 
  # scale_color_gradient2(low = "darkblue", mid = "white", high = "red", midpoint = 0, na.value = "white") + 
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "red", midpoint = 0, na.value = "white") + 
  coord_fixed(ratio = 6) +
  labs(x = "Gene", y = NULL) +
  # facet_grid(. ~ pathway, scales = "free") +
  theme_bw() + 
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
        # axis.text.x = element_blank(),
        panel.grid = element_blank())



