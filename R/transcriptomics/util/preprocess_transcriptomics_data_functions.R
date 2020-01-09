##-------------------------------------------------------------------------------------------------------##
# Methods for preprocessing data before further analysis

##-------------------------------------------------------------------------------------------------------##

# process data ahead of analysis
# remove probes with zero counts throughout

preprocess_chip_data <- function(count_data, gene_ids, org.Hs.eg.db, freq_thresh = 0.75)
{
nix <- which(rowSums(count_data) == 0) #find index of rows with all 0s
gene_ids_starting <- gene_ids
# remove gene ids that match with rows that contain all 0s
if(length(nix) != 0) {
  gene_ids <- gene_ids[-nix]
  count_data <- count_data[-nix, ]
}

# remove gene ids with multiple entrez mappings
gkey <- AnnotationDbi::select(x = org.Hs.eg.db, 
                              keys = gene_ids, 
                              columns = "ENTREZID", 
                              keytype = "SYMBOL")

# this basically looks for which genes appear more than once (i.e., have more than one entrez ID per gene ID)
nix <- which(sapply(gene_ids, function(y) length(which(gkey$SYMBOL == y))) != 1)
if(length(nix) != 0) {
  gene_ids <- gene_ids[-nix]
  count_data <- count_data[-nix, ]
}

# remove genes with NA entrez id
gkey <- AnnotationDbi::select(x = org.Hs.eg.db, 
                              keys = gene_ids, 
                              columns = "ENTREZID", 
                              keytype = "SYMBOL")

nix <- which(is.na(gkey$ENTREZID))
if(length(nix) != 0) {
  gene_ids <- gene_ids[-nix]
  count_data <- count_data[-nix, ]
}

# frequency across samples
# minimum 5 read counts, at least in 75% of the samples
x1 <- apply(count_data, 1, function(y) length(which(y >= 5)))
x1 <- which(x1 >= round(freq_thresh * ncol(count_data)))

count_data <- count_data[x1, ]
gene_ids <- gene_ids[x1]
genes_filtered<- setdiff(gene_ids_starting, gene_ids)

data_to_return <- list("count_data" = count_data, "gene_ids" = gene_ids, "genes_filtered" = genes_filtered)
return(data_to_return)
}

##-------------------------------------------------------------------------------------------------------##

# preprocess COPD external data
preprocess_copd_external_data <- function(count_data, gene_ids, org.Hs.eg.db)
{
  nix <- which(rowSums(count_data) == 0) #find index of rows with all 0s
  gene_ids_starting <- gene_ids
  # remove gene ids that match with rows that contain all 0s
  if(length(nix) != 0) {
    gene_ids <- gene_ids[-nix]
    count_data <- count_data[-nix, ]
  }
  
  # remove gene ids with multiple entrez mappings
  gkey_copd_all <- AnnotationDbi::select(org.Hs.eg.db,
                                         keys = gene_ids, 
                                         keytype = "ENSEMBL",
                                         columns=c("ENTREZID", "SYMBOL", "GENENAME")) # , "ENTREZID", "GENENAME"
  
  # ddply(df, "X1", numcolwise(sum))
  
  # this basically looks for which genes appear more than once (i.e., have more than one entrez ID per gene ID)
  nix <- which(sapply(gene_ids, function(y) length(which(gkey_copd_all$ENSEMBL == y))) != 1)
  if(length(nix) != 0) {
    gene_ids <- gene_ids[-nix]
    count_data <- count_data[-nix, ]
  }
  
  # remove genes with NA entrez id
  gkey_copd_semifilter <- AnnotationDbi::select(x = org.Hs.eg.db, 
                                                keys = gene_ids, 
                                                columns = c("ENTREZID", "SYMBOL", "GENENAME"), 
                                                keytype = "ENSEMBL")
  
  nix <- which(is.na(gkey_copd_semifilter$ENTREZID))
  if(length(nix) != 0) {
    gene_ids <- gene_ids[-nix]
    count_data <- count_data[-nix, ]
  }
  
  gkey_copd_semifilter2 <- AnnotationDbi::select(x = org.Hs.eg.db, 
                                                keys = gene_ids, 
                                                columns = c("ENTREZID", "SYMBOL", "GENENAME"), 
                                                keytype = "ENSEMBL")
  
  nix <- which(is.na(gkey_copd_semifilter2$SYMBOL))
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
  gene_ids_filtered<- setdiff(gene_ids_starting, gene_ids)
  gkey_copd_filtered <- AnnotationDbi::select(x = org.Hs.eg.db, 
                                                keys = gene_ids, 
                                                columns = c("ENTREZID", "SYMBOL", "GENENAME"), 
                                                keytype = "ENSEMBL")
  
  data_to_return <- list("count_data" = count_data, "gene_ids" = gene_ids, "genes_filtered" = gene_ids_filtered, "gkey_filtered" = gkey_copd_filtered)
  return(data_to_return)
}
##-------------------------------------------------------------------------------------------------------##

compute_avg_tissue_exp <- function(biogps_data)
{
  # average of tissues ----
  utis <- unique(biogps_data$tissues)
  avg_E <- matrix(0, nrow = nrow(biogps_data$expression), ncol = length(utis))
  for (i in seq(1, length(utis))) {
    u1 <- which(biogps_data$tissues == utis[i])
    if (length(u1) > 1) {
      avg_E[, i] <- rowMeans(biogps_data$expression[, u1])
    } else {
      avg_E[, i] <- biogps_data$expression[, u1]
    }
  }
  return(avg_E)
}


##-------------------------------------------------------------------------------------------------------##
