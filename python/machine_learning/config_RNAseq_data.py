# -*- coding: utf-8 -*-
"""
Created on Sun Jun 16 15:15:53 2019

@author: malcantar
"""
import pandas as pd
from scipy.stats import zscore


# load organ-on-chip RNA-seq data
# will need to update this once we get the next batch of data
# input: none
# output: raw count matrix, gene ids, sample metadata
def load_chip_data():

    # loading list that contains genes after the count data was preprocessed
    with open('../../data/gene_lists/gene_ids_chip_preprocessed.txt', 'r') as gene_ids:
        gene_ids = gene_ids.readlines()
    gene_ids_chip = [s.strip('\n') for s in gene_ids]

    # loading metadata and count data
    sample_metadata_chip = pd.read_csv('../../data/RNA-seq_data_chip/organ_on_chip_metadata.csv')
    count_data_chip = pd.read_csv('../../data/RNA-seq_data_chip/count_data_organ_on_chip_preprocessed.csv')

    # transposing matrix so features are columns and indices are sample names
    count_data_chip = count_data_chip.T
    count_data_chip.columns = gene_ids_chip

    return(count_data_chip, gene_ids_chip, sample_metadata_chip)

# apply z-score normalization to data
# input: raw count matrix
# output: z-score normalized count matrix
def standardize_data(count_data_df):
    count_data_zScore_df = count_data_df.apply(zscore)
    return count_data_zScore_df

# grab subset of original data such that we only have genes that were both differentially expressed and secreted / have signal peptides
# input: original count matrix
# output: dataframe with all samples, dataframe with only control_18 and virus_18 samples, genes that are secreted / differentially expressed
def grab_signal_secreted_genes(count_data_df, sample_metadata, time=18):

    if time==18:
        with open('../../data/gene_lists/virus_secreted_signal_chip.txt', 'r') as secreted_biomarker_genes:
            secreted_biomarker_genes = secreted_biomarker_genes.readlines()
        secreted_biomarker_genes = [s.strip('\n') for s in secreted_biomarker_genes]

        sample_to_group = dict(zip(sample_metadata.sample_name,
                                       sample_metadata.group_time))
        samples_control18_virus_18 = [key for key, value in sample_to_group.items() if 'control_18' in value.lower() or 'virus_18' in value.lower()]

        count_data_secreted_signal_df = count_data_df.copy()[secreted_biomarker_genes]

        count_data_ss_cont_vir_18_df = count_data_secreted_signal_df.copy().loc[samples_control18_virus_18]

        return(count_data_secreted_signal_df, count_data_ss_cont_vir_18_df, secreted_biomarker_genes)
    else:
        with open('../../data/gene_lists/virus_secreted_signal_chip_48hrs.txt', 'r') as secreted_biomarker_genes:
            secreted_biomarker_genes = secreted_biomarker_genes.readlines()
        secreted_biomarker_genes = [s.strip('\n') for s in secreted_biomarker_genes]

        sample_to_group = dict(zip(sample_metadata.sample_name,
                                       sample_metadata.group_time))
        samples_control48_virus_48 = [key for key, value in sample_to_group.items() if 'control_48' in value.lower() or 'virus_48' in value.lower()]

        count_data_secreted_signal_df = count_data_df.copy()[secreted_biomarker_genes]

        count_data_ss_cont_vir_48_df = count_data_secreted_signal_df.copy().loc[samples_control48_virus_48]

        return(count_data_secreted_signal_df, count_data_ss_cont_vir_48_df, secreted_biomarker_genes)
