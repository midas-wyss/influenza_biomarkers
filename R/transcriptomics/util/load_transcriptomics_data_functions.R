##-------------------------------------------------------------------------------------------------------##
# Loading data and reorganizing to make easier to handle during analysis

##------------------------------------------------------------------------------------------------------ ##

# This function will load RNA-seq count - and sample meta-data for organ-on-chip (obtained by Ingber lab) during the first RNA-seq analysis 
# Once the new dataset arrives, I should concatenate them first then run this pipeline

load_original_chip_data <- function()
{
  # load in count data and extract gene IDs
  count_data <- read_csv("../../../data/RNA-seq_data_chip/original_chip_data/STAR_Gene_Counts.csv")
  gene_ids <- count_data$Gene_ID
  count_data <- count_data %>% dplyr::select(., -Gene_ID)
  
  # load metadata and add sample IDs that match count matrix
  metadata <- readxl::read_xlsx("../../../data/RNA-seq_data_chip/original_chip_data/052918_RNAseq_Sample IDs.xlsx")
  newIDs <- gsub("5/2", "May2", gsub("5/1", "May1", gsub("-", "_", metadata$`Tube ID`))) 
  metadata <- metadata %>% 
    dplyr::mutate(., sample_name = newIDs) # adding a new column with new sample ID names
  
  # changing order of metadata so it matches count_data
  # essentially, match finds the indices that match the order of the first argument, 
  # we are then indexing metadata with these new indices such that the order matches
  metadata <- metadata[match(colnames(count_data), metadata$sample_name), ] 
  
  # shortened names for control and test groups
  group <- character(length = nrow(metadata))
  group[metadata$Sample == "Control"] <- "control"
  group[metadata$Sample == "Virus Treated"] <- "virus"
  group[metadata$Sample == "Poly IC Treated"] <- "poly_ic"
  
  # adding new column with group name (i.e., control, virus, poly_ic)
  metadata <- metadata %>% 
    dplyr::mutate(., group = group)
  
  # making new matrix with metadata
  sample_data <- tibble(sample_name = metadata$sample_name,
                        time_point = as.numeric(gsub("h", "", metadata$`Time point`)),
                        group = metadata$group)
  sample_data <- sample_data %>% 
    dplyr::mutate(., group_time = paste(group, time_point, sep = "_"))
  
  data_to_return <- list("count_data" = count_data, "metadata" = sample_data, "gene_ids" = gene_ids)
  
  return(data_to_return)
  
}

##------------------------------------------------------------------------------------------------------ ##
# load biogps data which displays gene expression in different tissue types

load_bioGPS <- function()
{
  load(file = "../../../data/2019-02-13_biogps.RData")
  return(biogps_data)
}

##------------------------------------------------------------------------------------------------------ ##

load_zaas <- function() {
  load("data/Zaas_etal_microarray/zaas_data_05152018.RData")
  data_to_return <- list("expression" = data[[1]], "genes" = data[[2]], "metadata" = data[[3]])
  return(data_to_return)
}

##------------------------------------------------------------------------------------------------------ ##
# load external COPD data

load_COPD_external <- function()
{
  COPD_count_data_original <- read_csv("../../../data/COPD_external_RNAseq/GSE124180_gene_count_table.csv")
  COPD_metadata <- read_csv("../../../data/COPD_external_RNAseq/COPD_metadata_external_final.csv")
  COPD_count_data_original <-COPD_count_data_original[-1,]
  COPD_gene_ids_EMSEMBL <- COPD_count_data_original$ENSEMBL_GENEID
  COPD_count_data_original <- COPD_count_data_original %>% dplyr::select(., -ENSEMBL_GENEID)
  COPD_metadata <- COPD_metadata[match(colnames(COPD_count_data_original), COPD_metadata$sample_name_original), ] 
  colnames(COPD_count_data_original) <- COPD_metadata$sample_name
  write.csv(COPD_count_data_original, file = "../../../data/COPD_external_RNAseq/copd_count_data_original_nogenes.csv")
  
  copd_count <- read_csv("../../../data/COPD_external_RNAseq/copd_count_data_original_nogenes.csv")[,-1]
  
  data_to_return <- list("copd_count_data"=copd_count, "genes"=COPD_gene_ids_EMSEMBL, "metadata" = COPD_metadata)
  return(data_to_return)
}

##------------------------------------------------------------------------------------------------------ ##
# new COPD | NHB RNAseq data

load_new_COPD_NHB_RNAseq_data<- function()
{
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









