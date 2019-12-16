# BENAM ET AL ----
# data on small airway exposed to cigarette smoke. 


# packages ----
library(pd.hg.u95av2, quietly = TRUE)
library(limma, quietly = TRUE)
library(tidyverse)
library(GEOquery)


# functions ----
source("/Volumes/HOME/scripts/r/expression_analysis/limma_dge.R")
source("/Volumes/HOME/scripts/r/wyss-compbio-functions.R")


# get data ----
load(file="data/public_data/2018-02-19_GSE87098_COPD-kambez.RData")

# data processing ----
# clean up probe sets
xx <- select(org.Hs.eg.db,
             keys = as.character(A_copd$GB_ACC),
             keytype = "REFSEQ",
             columns = c("ENTREZID","SYMBOL","GENENAME"))

nix <- which(xx$REFSEQ == "")
nix <- union(nix,which(is.na(xx$ENTREZID)))
E_copd <- E_copd[-nix, ]
A_copd <- A_copd[-nix, ]
gene_mappings_copd <- xx[-nix, ]
rm(xx,nix)

# summarize genes (keep unique entrez ids)
kp <- filter_genes(entrez_id = gene_mappings_copd$ENTREZID,
                   expression_data = E_copd)

E_copd <- E_copd[kp, ]
gene_mappings_copd <- gene_mappings_copd[kp, ]
A_copd <- A_copd[kp, ]

# differential expression ----
# compute differentially expressed genes
# comparing the 2 populations globally
healthy <- grep("healthy", S_copd$characteristics_ch1)
copd <- grep("COPD", S_copd$characteristics_ch1)

healthy_smoker <- intersect(grep("healthy", S_copd$characteristics_ch1), grep("cigarette", S_copd$characteristics_ch1.3))
healthy_nosmoker <- intersect(grep("healthy", S_copd$characteristics_ch1), grep("no smoke", S_copd$characteristics_ch1.3))

copd_smoker <- intersect(grep("COPD", S_copd$characteristics_ch1), grep("cigarette", S_copd$characteristics_ch1.3))
copd_nosmoker <- intersect(grep("COPD", S_copd$characteristics_ch1), grep("no smoke", S_copd$characteristics_ch1.3))

# limma bit
all_fits <- limma_dge(E_copd,
                      caseIds = list(healthy_smoker, copd_smoker, copd_nosmoker),
                      ctrIds = healthy_nosmoker) #<--- copd no smoker vs healthy no smoker

y1 <- limma_dge(E_copd,caseIds = copd_nosmoker,ctrIds = healthy_nosmoker) #<--- copd no smoker vs healthy no smoker
y2 <- limma_dge(E_copd,caseIds = copd_smoker,ctrIds = healthy_smoker) #<--- copd smoker vs healthy smoker
y3 <- limma_dge(E_copd,caseIds = copd_smoker,ctrIds = healthy_nosmoker) #<--- copd smoker vs healthy no smoker
y4 <- limma_dge(E_copd,caseIds = healthy_smoker,ctrIds = healthy_nosmoker) #<--- healthy smoker vs healthy no smoker
y5 <- limma_dge(E_copd,caseIds = copd_smoker,ctrIds = copd_nosmoker) #<--- copd smoker vs copd no smoker


# set of differentially expressed genes based on p-value threshold only
# diff_genes <- which(res$adj.P.Val < 0.01) # <-- less than a 1% false-discovery rate
z1 <- which(abs(y1$logFC) > 1 & y1$adj.P.Val < 0.01) # copd vs healthy, no smokers
z2 <- which(abs(y2$logFC) > 1 & y2$adj.P.Val < 0.01) # copd vs healthy, smokers
z3 <- which(abs(y3$logFC) > 1 & y3$adj.P.Val < 0.01) # copd smoker vs healhty non-smoker


# gene set enrichment with rpegeos ----
library(rpegeos)
copd_enr <- rpegeos::enrich_geneset(gene_mappings_copd$ENTREZID[z1])


dplyr::mutate(., probability_random = replace(probability_random, probability_random == 0, 1/10000)) %>% dplyr::mutate(., enr_score = -log10(probability_random) + 1 + (cosine_similarity * number_genes))

# DRUID COPD ----
copd_druid <- DRUID::concoct(dge_matrix = cbind(y1$logFC,y1$adj.P.Val), 
                             num_random = 10000, 
                             druid_direction = "neg", 
                             fold_thr = 1, 
                             pvalue_thr = 0.01, 
                             entrez = gene_mappings_copd$ENTREZID)



