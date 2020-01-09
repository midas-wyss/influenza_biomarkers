##------------------------------------------------------------------------------------------------------ ##
# Defining functions that will be used throughout analysis. In particular, these functions are used for 
# differential expression analysis and pathway enrichment

##------------------------------------------------------------------------------------------------------ ##
# Stouffer specificity 

# function for assessing tissue specificity
# INPUT: gene expression matrx (gene x sample) and
# OUTPUT: stouffer coefficient

stouffer_specificity <- function(expression_data) {
  # stouffer method:
  # (((x - x_hat) / sd(x)) + ((y - y_hat) / sd(y))) / sqrt(2)
  # where x is row and y is column
  
  means_rows <- rowMeans(expression_data) # average across genes
  means_cols <- colMeans(expression_data) # average across samples
  
  sd_rows <- apply(expression_data, 1, sd) # standard deviation of gene
  sd_cols <- apply(expression_data, 2, sd) # standard deviation of samples
  
  # calculate stoffer specificity
  z_rows <- apply(expression_data, 2, function(t) (t - means_rows) / sd_rows) # z-scores with respect to rows
  z_cols <- t(apply(expression_data, 1, function(t) (t - means_cols) / sd_cols)) # z-scores with respect to columns
  # z_cols <- t(z_cols)
  z_tot <- (z_rows + z_cols) / sqrt(2)
  return(z_tot) 
}

##------------------------------------------------------------------------------------------------------ ##
# Differential expression using limma -- used on microarray data

# performs differential expression on microarray expression data
# INPUT: gene expression matrix (gene x sample), indices for cases, indices for controls
# OUTPUT: fold-change matrix (with associated p-values)

# function for differential expression with limma
# takes in expression matrix, case IDs, and control IDs
limma_dge <- function(expression_data,caseIds=list(),ctrIds)
{
  require(limma)
  if(!is.list(caseIds)) #ie, only one case
  {
    tmp_data <- expression_data[,c(ctrIds,caseIds)]
    
    expdes <- matrix(0,nrow=ncol(tmp_data),2) # experimental design -- initialize matrix with 0s -- these will eventually become one-hot encodings for case vs. control
    expdes[seq(1,length(ctrIds)),1] <- 1
    expdes[seq(length(ctrIds)+1,length(ctrIds)+length(caseIds)),2] <- 1
    colnames(expdes) <- c("control","case")
    
    contmat <- makeContrasts(case-control,levels=expdes) # contrast matrix
    
    limma.fit <- lmFit(tmp_data,design=expdes)
    limma.fit <- contrasts.fit(limma.fit,contmat)
    limma.fit <- eBayes(limma.fit)
    
    res <- topTable(limma.fit,1,number=nrow(tmp_data),sort.by = "none",adjust.method = "fdr")
    return(res)
    
  } else {
    # one-hot encode
    expdes <- matrix(0,nrow=ncol(expression_data),1) 
    expdes[ctrIds,1] <- 1 # controls
    
    for(i in seq(1,length(caseIds)))
    {
      x <- matrix(0,nrow=ncol(expression_data),1) 
      x[caseIds[[i]],1] <- 1
      expdes <- cbind(expdes,x)
    }
    
    tmp_data <- expression_data[,c(ctrIds,unlist(caseIds))]
    expdes <- expdes[c(ctrIds,unlist(caseIds)),]
    
    colnames(expdes) <- c("control",paste("case",seq(1,length(caseIds)),sep=""))
    myContrasts <- c(paste(colnames(expdes)[2:ncol(expdes)],"-",colnames(expdes)[1],sep=""))
    contmat <- contmat <- eval(as.call(c(as.symbol("makeContrasts"),as.list
                                         (myContrasts),levels=list(expdes)))) # contrast matrix
    
    limma.fit <- lmFit(tmp_data,design=expdes)
    limma.fit <- contrasts.fit(limma.fit,contmat)
    limma.fit <- eBayes(limma.fit)
    
    return(limma.fit)
  }
}

##------------------------------------------------------------------------------------------------------ ##
## Gene set enrichment analysis using the GAGE package (https://bioconductor.org/packages/release/bioc/html/gage.html)

# takes in matrix with fold changes and gene names, and returns pathway enrichment figures and list of up- and down-regulated genes
GSEA_Gage <- function(res_all, DEmethod, G = NULL, output_loc, time_point, draw_map = FALSE) {
  if(DEmethod == "DESeq2")
  {
    foldchanges = res_all$log2FoldChange
    names(foldchanges) = res_all$entrez
    # head(foldchanges)
  }
  
  else if(DEmethod == "limma") 
  {
    foldchanges = res_all$logFC
    names(foldchanges) = G$ENTREZID
  }
  
  else
  {
    stop(cat(paste("Invalid DE method. Please enter either", shQuote("DESeq2"), "or",shQuote("limma"), sep="")))
  }
  
  keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)
  lapply(keggres, head)
  
  keggrespathways_up = data.frame(id=rownames(keggres$greater), keggres$greater) %>% 
    tbl_df() %>% 
    filter(row_number()<=20) %>% 
    .$id %>% 
    as.character()
  
  keggrespathways_maps_up = data.frame(id=rownames(keggres$greater), keggres$greater) %>% 
    tbl_df() %>% 
    filter(row_number()<=5) %>% 
    .$id %>% 
    as.character()

  # Get the pathways
  keggresids_up = substr(keggrespathways_maps_up, start=1, stop=8)
  
  if(draw_map){
  # Define plotting function for applying later
  plot_pathway_up = function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa", new.signature=FALSE, kegg.dir = "../figs/pathway_maps/")
  
  # plot multiple pathways (plots saved to disk and returns a throwaway list object)
  tmp = sapply(keggresids_up, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa", , kegg.dir = "../figs/pathway_maps/"))
}
  keggrespathways_down = data.frame(id=rownames(keggres$less), keggres$less) %>% 
    tbl_df() %>% 
    filter(row_number()<=20) %>% 
    .$id %>% 
    as.character()
  keggresids_down = substr(keggrespathways_down, start=1, stop=8)
  
  upregulated <- subset(keggres$greater, rownames(keggres$greater) %in% keggrespathways_up)
  downregulated <- subset(keggres$less, rownames(keggres$less) %in% keggrespathways_down)
  
  
  print(upregulated %>%tbl_df%>% ggplot(aes(y = stat.mean , x = forcats::fct_reorder(substr(rownames(upregulated),start=10, stop=10e2), stat.mean), size = set.size, colour= p.val)) +
    geom_point()+
    labs(x = "Gene set", y = "Enrichment score", colour="p value", size="Count") +   coord_flip() +  theme_bw() + 
    theme(axis.title = element_text(size = 16, color = "black"),
          axis.text = element_text(size = 16, color = "black") ))
  
  ggsave(filename = paste(output_loc,"_",time_point,"_upregulated_GAGE.pdf",sep=""),
         plot = last_plot(),
         device = "pdf",
         height = 8,
         width = 10,
         scale = 1,
         dpi = 400
  )

  ggsave(filename =paste(output_loc,"_",time_point,"_upregulated_GAGE.png",sep=""),
         plot = last_plot(),
         device = "png",
         height = 8,
         width = 10,
         scale = 1,
         dpi = 400
  )
  
  # ggsave(filename = paste(output_loc,"_upregulated_GAGE.svg",sep=""), 
  #        plot = last_plot(), 
  #        device = "svg", 
  #        height = 8,
  #        width = 10,
  #        scale = 1,
  #        dpi = 400
  # )
  
  print(downregulated %>%tbl_df%>% ggplot(aes(y = stat.mean , 
                                              x = forcats::fct_reorder(substr(rownames(downregulated),start=10, 
                                              stop=10e2), stat.mean), size = set.size, colour= p.val)) + geom_point()+
    labs(x = "Gene set", y = "Enrichment score", colour="p value", size="Count") +   coord_flip() +  theme_bw() + 
    theme(axis.title = element_text(size = 16, color = "black"),
          axis.text = element_text(size = 16, color = "black") ))
  
  # ggsave(filename = paste(output_loc,"_downregulated_GAGE.svg",sep=""), 
  #        plot = last_plot(), 
  #        device = "svg", 
  #        height = 8,
  #        width = 10,
  #        scale = 1, 
  #        dpi = 400 
  # )
  # 
  ggsave(filename = paste(output_loc,"_",time_point,"_downregulated_GAGE.png",sep=""), 
         plot = last_plot(), 
         device = "png", 
         height = 8,
         width = 10,
         scale = 1, 
         dpi = 400 
  )
  
  ggsave(filename = paste(output_loc,"_",time_point, "_downregulated_GAGE.pdf",sep=""), 
         plot = last_plot(), 
         device = "pdf", 
         height = 8,
         width = 10,
         scale = 1, 
         dpi = 400 
  )
  
  result <- list("upregulated" = upregulated, "downregulated" = downregulated) 
  return(result) 
  
}

##------------------------------------------------------------------------------------------------------ ##

# perform differential expression analysis on specified groups
compare_group_DESeq2 <- function(dds, contrast_groups, gene_ids, MHT = "fdr")
  {
  res <- results(dds,
                     contrast = contrast_groups,
                     pAdjustMethod = MHT,
                     cooksCutoff = FALSE)
  
  res <- res %>% 
    tibble::add_column(., gene = gene_ids, .before = 1)
  
  return(res)
}

##------------------------------------------------------------------------------------------------------ ##
# gene set enrichment with rpegeos ----
# gnames <- gene_ids[virus_specific]

rpegeos_enrichment <- function(genes, res, num_paths=20, method = "DESeq2", org.Hs.eg.db = org.Hs.eg.db, output_loc)
{
gkey <- AnnotationDbi::select(x = org.Hs.eg.db, 
                              keys = genes, 
                              columns = "ENTREZID", 
                              keytype = "SYMBOL")

if(method == "DESeq2")
{
  res_specific <- res[res$gene %in% genes,]
  chip_gset <- cbind(gkey$ENTREZID, res_specific$log2FoldChange, res_specific$padj)
  chip_gset <- apply(chip_gset, 2, as.numeric)
  chip_enr <- rpegeos::enrich_geneset(gene_set = chip_gset)
}
else
{
  res_specific <- res[res$gene_symbol  %in% genes,]
  chip_gset <- cbind(gkey$ENTREZID, res_specific$fold_change , res_specific$fdr_pvalue)
  chip_gset <- apply(chip_gset, 2, as.numeric)
  chip_enr <- rpegeos::enrich_geneset(gene_set = chip_gset)
}
print(chip_enr %>%
  dplyr::slice(1:num_paths) %>%
  ggplot() +
  geom_point(aes(y = enrichment_score, x = forcats::fct_reorder(geneset, enrichment_score), size = number_genes, color = probability_random), alpha = 1.0
  ) +
  labs(x = "Gene set", y = "Enrichment score", colour= "p value", size="Count") +
  coord_flip() + 
  theme_bw() + 
  theme(axis.title = element_text(size = 16, color = "black"),
        axis.text = element_text(size = 16, color = "black")))

ggsave(filename = paste(output_loc,".pdf",sep=""), 
       plot = last_plot(), 
       device = "pdf", 
       height = 8,
       width = 10,
       scale = 1, 
       dpi = 400)

ggsave(filename = paste(output_loc,".png",sep=""), 
       plot = last_plot(), 
       device = "png", 
       height = 8,
       width = 10,
       scale = 1, 
       dpi = 400)

# ggsave(filename = paste(output_loc,".svg",sep=""), 
#        plot = last_plot(), 
#        device = "svg", 
#        height = 8,
#        width = 10,
#        scale = 1, 
#        dpi = 400)

return(chip_enr)
}



##------------------------------------------------------------------------------------------------------ ##
# lung-specific flu biomarkers

lung_specific_biomarkers <- function(avg_tissue_exp,virus_specific_genes, biogps_data, stouf_thr=3, num_tissues=2, output_dir)
{
  ts_stouffer <- stouffer_specificity(avg_tissue_exp) # <-- see how i calculate the average of tissue expression
  print(hist(ts_stouffer, breaks = 50))
  ggsave(filename = paste(output_dir,"stouffer_histogram_biogps",".pdf",sep=""), 
         plot = last_plot(), 
         device = "pdf", 
         height = 8,
         width = 10,
         scale = 1, 
         dpi = 400)
  
  ggsave(filename = paste(output_dir,"stouffer_histogram_biogps",".png",sep=""), 
         plot = last_plot(), 
         device = "png", 
         height = 8,
         width = 10,
         scale = 1, 
         dpi =400)
  # ggsave(filename = paste("../../../figs/tissue_specificty/stouffer_histogram_biogps",".svg",sep=""), 
  #        plot = last_plot(), 
  #        device = "svg", 
  #        height = 8,
  #        width = 10,
  #        scale = 1, 
  #        dpi =400)
  # 
  
  tissue_spec <- tibble::tibble(count = c(length(which(apply(ts_stouffer, 1, function(y) length(which(abs(y) > 1))) != 0)),
                                        length(which(apply(ts_stouffer, 1, function(y) length(which(abs(y) > 2))) != 0)),
                                        length(which(apply(ts_stouffer, 1, function(y) length(which(abs(y) > 3))) != 0)),
                                        length(which(apply(ts_stouffer, 1, function(y) length(which(abs(y) > 4))) != 0)),
                                        length(which(apply(ts_stouffer, 1, function(y) length(which(abs(y) > 5))) != 0)),
                                        length(which(apply(ts_stouffer, 1, function(y) length(which(abs(y) > 6))) != 0))),
                              group = c("z > 1", "z > 2", "z > 3", "z > 4", "z > 5", "z > 6"))
  
  
  print(tissue_spec %>% ggplot() + geom_col(aes(y = count, x = group)))
  
  ggsave(filename = paste(output_dir, "tissue_stouffer_score_biogps",".pdf",sep=""), 
         plot = last_plot(), 
         device = "pdf", 
         height = 8,
         width = 10,
         scale = 1, 
         dpi = 400)
  
  ggsave(filename = paste(output_dir, "tissue_stouffer_score_biogps",".png",sep=""), 
         plot = last_plot(), 
         device = "png", 
         height = 8,
         width = 10,
         scale = 1, 
         dpi = 400)
  
  # ggsave(filename = paste("../figs/tissue_stouffer_score_biogps",".svg",sep=""), 
  #        plot = last_plot(), 
  #        device = "pdf", 
  #        height = 8,
  #        width = 10,
  #        scale = 1, 
  #        dpi = 400)
  
  stouf_thr <- stouf_thr
  num_tissues <- num_tissues
  
  spec_mat <- matrix(0, nrow = nrow(x = ts_stouffer), ncol = ncol(ts_stouffer))
  spec_mat[ts_stouffer >= stouf_thr] <- 1
  
  spec_genes <- which(rowSums(spec_mat) <= num_tissues)
  lung_genes <- biogps_data$genes$SYMBOL[intersect(which(spec_mat[, 2] == 1), spec_genes)]
  sort(lung_genes)
  # length(lung_genes)
  
  # intersect of viral biomarkers and lung specific 
  lung_specific_biomarkers_list <- list("lung_spcific" = intersect(lung_genes,virus_specific_genes), "lung_genes" = lung_genes)
  # print(lung_specific_biomarkers)
  return(lung_specific_biomarkers_list)
}

##------------------------------------------------------------------------------------------------------ ##
# differential expression on Zaas data

DE_limma_Zaas <- function(S, E)
{
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
  
  results_to_return <- list("flu" = xx1, "rsv" = xx2, "hrv" = xx3)
  
  return(results_to_return)
  
}
  
 
##------------------------------------------------------------------------------------------------------ ##
# network stuff
network_analysis <- function(DN,ttrust_data, org.Hs.eg.db, gene_ids, virus_specific, lung_genes)
{
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
    dplyr::filter(., x %in% xx$ENTREZID | y %in% xx$ENTREZID) %>% 
    dplyr::filter(., y %in% xx$ENTREZID | x %in% xx$ENTREZID) %>% 
    dplyr::filter(., change_type != "edge present, same sign") %>% 
    dplyr::mutate(., x_lung = x_name %in% lung_genes) %>%
    dplyr::mutate(., x_lung = replace(x_lung, x_lung == TRUE, 1), x_lung = replace(x_lung, x_lung == FALSE, 0)) %>%
    dplyr::mutate(., y_lung = y_name %in% lung_genes) %>%
    dplyr::mutate(., y_lung = replace(y_lung, y_lung == TRUE, 1), y_lung = replace(y_lung, y_lung == FALSE, 0)) %>% 
    dplyr::filter(., x_lung == 1 | y_lung == 1) %>% 
    dplyr::filter(., x_name %in% gene_ids, y_name %in% gene_ids)
  
  
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
  data_to_return <- list("node_names" = node_names, "tmp" = tmp)
  return(data_to_return)
}


##------------------------------------------------------------------------------------------------------ ##
## Compare paths GAGE (Zaas and chip)
compare_chip_Zaas <- function(enriched_paths_chip, enriched_paths_flu, direction = "upregulated", output_dir)
  {
  if(direction == "upregulated")
  {
    enriched_paths_cp <-enriched_paths_chip$upregulated 
    enriched_paths_fl <- enriched_paths_flu$upregulated
  }
  else
  {
    enriched_paths_cp <-enriched_paths_chip$downregulated 
    enriched_paths_fl <- enriched_paths_flu$downregulated
  }
  g1 <- enriched_paths_cp %>%tbl_df %>% dplyr::slice(1:20)
  g2 <- enriched_paths_fl %>%tbl_df%>% dplyr::slice(1:20)

  g1 <- g1 %>% 
    tibble::add_column(., group = "organ on chip") %>% 
    tibble::add_column(., rank = seq(1, nrow(g1))) %>%
    dplyr::mutate(., norm_score = (stat.mean - min(stat.mean)) / (max(stat.mean) - min(stat.mean)))
  
  g2 <- g2 %>% 
    tibble::add_column(., group = "zaas et al") %>% 
    tibble::add_column(., rank = seq(1, nrow(g2))) %>%
    dplyr::mutate(., norm_score = (stat.mean - min(stat.mean)) / (max(stat.mean) - min(stat.mean)))
  
  g4 <- rbind(g1, g2)
  g4 <- add_column(g4, pathways = c(substr(rownames(enriched_paths_cp), 
                                           start=10, stop=10000),
                                    substr(rownames(enriched_paths_fl), 
                                           start=10, stop=10000)))
  print(g4 %>% 
    #dplyr::filter(., group != "network data") %>%
    ggplot() + 
    geom_point(aes(x = forcats::fct_reorder(pathways, stat.mean), y = stat.mean, color = group, size = set.size)) +
    scale_color_manual(values = c("#ff9933", "#cc99cc")) + # , "#99cc99"
    labs(x = "Gene set", y = "Enrichment score") +
    coord_flip() + 
    facet_grid(. ~ group, scales = "free") +
    theme_bw() + 
    theme(axis.text = element_text(size = 14,color = "black"),
          axis.title = element_text(size = 14,color = "black"),
          strip.text = element_text(size = 14,color = "black"),
          legend.position = "none"))
  
  ggsave(filename = paste(output_dir, "GSEA_GAGE_comparison_chip_zaas","_",direction,".pdf",sep=""), 
         plot = last_plot(), 
         device = "pdf", 
         height = 8,
         width = 10,
         scale = 1, 
         dpi = 400)
  ggsave(filename = paste(output_dir, "GSEA_GAGE_comparison_chip_zaas","_",direction,".png",sep=""), 
         plot = last_plot(), 
         device = "png", 
         height = 8,
         width = 10,
         scale = 1, 
         dpi = 400)
  
  # ggsave(filename = paste("../figs/GSEA_GAGE_comparison_chip_zaas","_",direction,".svg",sep=""), 
  #        plot = last_plot(), 
  #        device = "svg", 
  #        height = 8,
  #        width = 10,
  #        scale = 1, 
  #        dpi = 400)
}


##------------------------------------------------------------------------------------------------------ ##
## compare chip and zaas data from rpegeos

compare_zaas_chip_rpegeos <- function(flu_enr, chip_enr, numPaths,output_dir)
{
a1 <- flu_enr %>% dplyr::slice(1:numPaths)
#a2 <- subnet_enr %>% dplyr::slice(1:numPaths)
a3 <- chip_enr %>% dplyr::slice(1:numPaths)

a1 <- a1 %>% 
  tibble::add_column(., group = "zaas et al. data") %>% 
  tibble::add_column(., rank = seq(1, nrow(a1))) %>%
  dplyr::mutate(., norm_score = (enrichment_score - min(enrichment_score)) / (max(enrichment_score) - min(enrichment_score)))

#a2 <- a2 %>% 
# tibble::add_column(., group = "network data") %>% 
#tibble::add_column(., rank = seq(1, nrow(a2))) %>%
# dplyr::mutate(., norm_score = (enrichment_score - min(enrichment_score)) / (max(enrichment_score) - min(enrichment_score)))

a3 <- a3 %>% 
  tibble::add_column(., group = "organ-on-chip data") %>% 
  tibble::add_column(., rank = seq(1, nrow(a3))) %>%
  dplyr::mutate(., norm_score = (enrichment_score - min(enrichment_score)) / (max(enrichment_score) - min(enrichment_score)))

a4 <- rbind(a1, a3) # a2

print(a4 %>% 
  #dplyr::filter(., group != "network data") %>%
  ggplot() + 
  geom_point(aes(x = forcats::fct_reorder(geneset, enrichment_score), y = enrichment_score, color = group, size = number_genes)) +
  scale_color_manual(values = c("#ff9933", "#cc99cc")) + # , "#99cc99"
  labs(x = "Gene set", y = "Enrichment score") +
  coord_flip() + 
  facet_grid(. ~ group, scales = "free") +
  theme_bw() + 
  theme(axis.text = element_text(size = 14,color = "black"),
        axis.title = element_text(size = 14,color = "black"),
        strip.text = element_text(size = 14,color = "black"),
        legend.position = "none"))

ggsave(filename = paste(output_dir, "GSEA_rpegeos_compare_chip_zaas.pdf", sep=""), 
       plot = last_plot(), 
       device = "pdf", 
       height = 8,
       width = 10,
       scale = 1, 
       dpi = 400
)

ggsave(filename = paste(output_dir, "GSEA_rpegeos_compare_chip_zaas.png", sep=""), 
       plot = last_plot(), 
       device = "png", 
       height = 8,
       width = 10,
       scale = 1, 
       dpi = 400
)

# ggsave(filename = "../figs/GSEA_rpegeos_compare_chip_zaas.svg", 
#        plot = last_plot(), 
#        device = "svg", 
#        height = 8,
#        width = 10,
#        scale = 1, 
#        dpi = 400 
# )
}

##------------------------------------------------------------------------------------------------------ ##
## DE for Huang et al. data
DE_Huang_Woods <- function(S, E, reverse = FALSE)
  {
  if(!reverse){
  symp <- which(list(S$disease)[[1]] == "Symptomatic")
  asymp <- which(list(S$disease)[[1]] == "Asymptomatic")
  }
  else {
    asymp <- which(list(S$disease)[[1]] == "Symptomatic")
    symp <- which(list(S$disease)[[1]] == "Asymptomatic")
  }
  xx1 <- limma_dge(expression_data = E,
                   ctrIds = asymp,
                   caseIds = symp)
  return(xx1)
  
  }
  
##------------------------------------------------------------------------------------------------------ ##
# Create volcano plot

plot_volcano <- function(res, x_min, x_max, pcutoff = 0.05, output_name)
{
  EnhancedVolcano(res,
                  lab = res$gene, # res_c48_c18
                  x = 'log2FoldChange',
                  y = 'padj',
                  pCutoff = pcutoff,
                  FCcutoff = 1.0,
                  xlim = c(-2, 7))
  ggsave(filename = paste(output_name, '.pdf',sep=""),
         plot = last_plot(), 
         device = "pdf", 
         height = 8,
         width = 10,
         scale = 1, 
         dpi = 400 
  )
  ggsave(filename = paste(output_name,'.png',sep=""),
         plot = last_plot(), 
         device = "png", 
         height = 8,
         width = 10,
         scale = 1, 
         dpi = 400 
  )
  
} 
  
  