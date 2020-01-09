# -*- coding: utf-8 -*-
"""
Created on Sun Jun 16 15:22:00 2019

@author: malcantar
"""
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from sklearn import decomposition
from pylab import *


# run PCA on normalized count matrix
# input: normalized count matrix
# output: PCA matrix and PCA object

def run_PCA(count_data_zScore_df):

    # creating list of principal component names to simplify column naming
    column_names = []
    for numPC in range(1, min(count_data_zScore_df.shape) + 1):
            column_names.append('principal component ' + str(numPC))

    # computing principal components
    pca_obj = decomposition.PCA(n_components=numPC)
    pca_mat= pca_obj.fit_transform(count_data_zScore_df)
    pca_df = pd.DataFrame(data = pca_mat,
    columns = column_names)
    pca_df.index = count_data_zScore_df.index # index names will be sample names
    return(pca_df, pca_obj)

# computes PCA loadings
# inputs: PCA object and gene ids (needed to map PCA index to gene)
# output: loading matrix
def compute_pca_loadings(pca_obj, gene_ids):

    # pca loadings calculation
    pca_loadings = pca_obj.components_.T * np.sqrt(pca_obj.explained_variance_) # calculating loadings

    # again, creating list of column names to simplify naming columns
    column_names = []
    for num_load in range(1, pca_loadings.shape[1] + 1):
        column_names.append('principal component ' + str(num_load))

    # placing loadings in a dataframe
    loadings_PCA_df = pd.DataFrame(data = pca_loadings,
                              columns = column_names)
    # making gene IDs the indices
    loadings_PCA_df.index = gene_ids
    return(loadings_PCA_df)

 # create PCA plot
# inputs: PCA matrix, PCA object, metadata, gene IDS, PCs you want to plot
# output: PCA plot with principal components specificed in argument 'PCs'
def plotPCA(pca_df, pca_obj, sample_metadata, gene_ids_chip,
            PCs, save_plot=False, output_loc=False):

    #####--------------------------------adding groups to count data matrix--------------------------------#####

    sample_to_group = dict(zip(sample_metadata.sample_name,
                               sample_metadata.group_time))
    for sample in pca_df.index:
         pca_df.loc[sample,'group'] = sample_to_group[sample]

    #####-------------------------------getting PCs you want to plot--------------------------------#####

    PCa = PCs[0]
    PCb = PCs[1]

    #####-------------------------------plotting PCs--------------------------------#####

    # plotting both specified PCs colored by presence-absence of top genes
    varianceExplainedAbs = pca_obj.explained_variance_ratio_ * 100
    fig = plt.figure(figsize = (8,8))
    ax = fig.add_subplot(1,1,1)
    ax.set_xlabel('Principal Component %s (%0.2f)' % (str(PCa), varianceExplainedAbs[PCa-1]), fontsize = 20)
    ax.set_ylabel('Principal Component %s (%0.2f)' % (str(PCb), varianceExplainedAbs[PCb-1]), fontsize = 20)
    ax.set_title('2 component PCA', fontsize = 24)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)

    # conditions extracted based on the groups beings tested
    conditions = set([value for key, value in sample_to_group.items() if key in list(pca_df.index)])

    # looping through and coloring each allele combination a different color
    colors = sns.color_palette()[0:len(list(set(sample_metadata.group_time)))]

    for target, color in zip(conditions, colors):
        indicesToKeep = pca_df['group'] == target
        ax.scatter(pca_df.loc[indicesToKeep, 'principal component ' + str(PCa)]
                   , pca_df.loc[indicesToKeep, 'principal component '  + str(PCb)]
                    , c = color
                    , s = 50)
    ax.legend(conditions, prop={'size': 18})

    if save_plot and output_loc:
        plt.savefig(output_loc + 'PC'+ str(PCa) + '_' + 'PC' + str(PCb) + '.pdf', dpi=400)
        plt.savefig(output_loc + 'PC'+ str(PCa) + '_' + 'PC' + str(PCb) + '.png', dpi=400)
        plt.savefig(output_loc + 'PC'+ str(PCa) + '_' + 'PC' + str(PCb) + '.svg', dpi=400)
# plots PCA loadings for specified PCs
# input: PCA object, list of 2 loadings you want to plot
# output: plot of loadings for whataver PCs you want to analyze
def plot_loadings(pca_obj, secreted_biomarkers, loadings_to_plot):

    # computing loadings using previously defined function
    loadings = compute_pca_loadings(pca_obj, secreted_biomarkers)

    # plotting loadings in order of descending value
    figure(num = 1, figsize=(15, 6), dpi = 400)
    plt.subplot(1, 2, 1)
    loadings_sorted = loadings.sort_values('principal component 1',ascending=False)
    ax = loadings_sorted['principal component 1'].plot(kind='bar',facecolor='#AA0000')
    # ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
    ax.patch.set_facecolor('#FFFFFF')
    ax.spines['bottom'].set_color('#CCCCCC')
    ax.spines['bottom'].set_linewidth(1)
    ax.spines['left'].set_color('#CCCCCC')
    ax.spines['left'].set_linewidth(1)
    plt.title('PC1 loadings', fontsize=24)

    plt.subplot(1, 2, 2)
    loadings_sorted = loadings.sort_values('principal component 2',ascending=False)
    ax = loadings_sorted['principal component 2'].plot(kind='bar',facecolor='blue')
    # ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
    ax.patch.set_facecolor('#FFFFFF')
    ax.spines['bottom'].set_color('#CCCCCC')
    ax.spines['bottom'].set_linewidth(1)
    ax.spines['left'].set_color('#CCCCCC')
    ax.spines['left'].set_linewidth(1)
    plt.title('PC2 loadings', fontsize=24)
    plt.show()
