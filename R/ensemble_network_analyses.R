load("data/2019-05-09_biogps_network_consensus.RData")

# network analyses of virus specific genes ----
vs <- gkey[virus_specific, ]

vs_net <- CNET$consensus_network %>% 
  dplyr::filter(., x %in% vs$SYMBOL | y %in% vs$SYMBOL) %>% 
  dplyr::filter(., absolute_majority == 1)


# network analyses of lung specific genes ----
lspec <- lung_genes

lspec_net <- CNET$consensus_network %>% 
  dplyr::filter(., x %in% lspec | y %in% lspec) %>% 
  dplyr::filter(., super_majority == 1)


# looking at ttrust data
a1 <- nsmblR::ttrust_data[which(nsmblR::ttrust_data$tf %in% vs$SYMBOL), ]
a2 <- nsmblR::ttrust_data[which(nsmblR::ttrust_data$target %in% vs$SYMBOL), ]

intersect(a1$target, a2$target)

b1 <- nsmblR::ttrust_data[which(nsmblR::ttrust_data$tf %in% lspec), ]
b2 <- nsmblR::ttrust_data[which(nsmblR::ttrust_data$target %in% lspec), ]
