##------------------------------------------------------------------------------------------------------ ##
# new COPD | NHB RNAseq data

load_data<- function() {
  # load raw metadata
  metadata_raw_copd_nhb <- readxl::read_xlsx("data/RNA-seq_data_chip/COPD_NHB_data/RNAseq_NHB_COPD_metadata_new.xlsx")
  
  # function that will help create new column names with form condition#_time_donorID
  reorder_names <- function(x)
  {
    x_split <- unlist(strsplit(x, " "))
    x_len <- length(x_split)
    
    if(length(x_split[1:(x_len-2)]) > 1)
    {
      return(paste(x_split[x_len], x_split[(x_len-1)], x_split[1], x_split[2], sep="_"))
    }
    else
    {
      return(paste(x_split[x_len], x_split[(x_len-1)], x_split[1], sep="_"))
    }
  }
  
  # create final metadata csv file
  new_sample_names <-unlist(lapply(metadata_raw_copd_nhb$sample_name, reorder_names))
  metadata_copd_nhb <- metadata_raw_copd_nhb %>% 
    dplyr::mutate(., new_sample_names = new_sample_names,
                  group_time = paste(group, time, sep = "_"),
                  condition_time = paste(condition, time, sep = "_"),
                  condition_group = paste(condition, group, sep = "_"),
                  condition_group_time = paste(condition, group, time, sep="_"))

  # uncomment if you want to save new csv file
  # write_csv(metadata_copd_nhb, "../../../data/RNA-seq_data_chip/COPD_NHB_data/metadata_copd_nhb.csv")
  
  # load in count data and extract gene IDs
  count_data <- read_csv("data/RNA-seq_data_chip/COPD_NHB_data/RP7193_hg19_STAR_raw_gene_counts.csv")
  gene_ids <- count_data$Gene_ID
  count_data <- count_data %>% dplyr::select(., -Gene_ID)
  
  # drop well that did not run in RNA seq pipeline
  count_data <- count_data[,!names(count_data) %in% 's2-G3']
  metadata_copd_nhb <- metadata_copd_nhb[!metadata_copd_nhb$well_name %in% 's2-G3',]
  
  # changing order of metadata so it matches count_data
  # essentially, match finds the indices that match the order of the first argument, 
  # we are then indexing metadata with these new indices such that the order matches
  metadata_copd_nhb <- metadata_copd_nhb[match(colnames(count_data), metadata_copd_nhb$well_name), ]

  # conditions we want to compare
  data_to_return <- list("count_data" = count_data, "metadata" = metadata_copd_nhb, "gene_ids" = gene_ids)
  
  return(data_to_return)
  
}


preprocess_data <- function(count_data, gene_ids, org.Hs.eg.db, freq_thresh = 0.75) {
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


##------------------------------------------------------------------------------------------------------ ##
# load biogps data which displays gene expression in different tissue types

load_bioGPS <- function()
{
  load(file = "data/2019-02-13_biogps.RData")
  return(biogps_data)
}









