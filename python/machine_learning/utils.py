#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 19:44:21 2019

@author: malcantar
"""

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from pylab import *


def summary_stats(ens):
    ens = ens.dropna(axis=0)
    cols_to_med = ['roc_auc', 'AP', 'accuracy'] 
    
    summary = pd.DataFrame()
    
    for col in cols_to_med:
        vals = ens[col].values
        med_idx = np.argsort(vals)[len(vals)//2]
        median = vals[med_idx]
        ci = np.percentile(vals, (2.5, 97.5))
        
        if np.isnan(median): #if nan
            continue
                
        summary[col+'_ci'] = [ci]
        summary[col+'_median'] = median
        
        
        med_row = ens.iloc[med_idx, :]
        if (col == 'roc_auc'):
            fpr = med_row['fpr']
            tpr = med_row['tpr']
            temp = pd.DataFrame({'fpr': [fpr], 'tpr': [tpr]})
            summary = pd.merge(summary, temp, how='left', 
                               left_index=True, right_index=True)
        if (col == 'AP'):
            precision = med_row['precision']
            recall = med_row['recall']
            temp = pd.DataFrame({'precision': [precision],
                                 'recall': [recall]})
            summary = pd.merge(summary, temp, how='left', 
                               left_index=True, right_index=True)
    return summary
        
def get_weights(weights_ens):
    weights_ens = weights_ens.abs()
    weights_ens = weights_ens.div(weights_ens.sum(axis=1), axis=0)
    weights_ens = weights_ens.sum(axis=0).transpose()
    weights_ens = pd.DataFrame(data=weights_ens, 
                               columns=['Feature importance'])
    #weights_ens = weights_to_pathways(weights_ens, chemfile)
    weights_ens = weights_ens.sort_values(by='Feature importance',
                                          ascending=False)
    return weights_ens
         
def to_array(str_array): 
    float_array = [float(string) for string in str_array[1:-1].split()]
    return float_array

# one-hot encode control (0) vs. virally-infected (1)
def one_hot_encode(count_data_zScore_df, sample_metadata_df):
    sample_to_group = dict(zip(sample_metadata_df.sample_name, 
                                       sample_metadata_df.group_time))

    conditions = [value for key, value in sample_to_group.items() if key in list(count_data_zScore_df.index)]
    virus_control_labels = [0 if 'control' in label else 1 for label in conditions]
    
    return(virus_control_labels)
    
# plot bar plot of feature weights
def plot_weights(weights_df):
    # plotting loadings in order of descending value
    figure(num = 1, figsize=(15, 6), dpi = 400)  
    plt.subplot(1, 2, 1)
    weights_sorted_df = weights_df.sort_values('Feature importance',ascending=False)
    ax = weights_sorted_df['Feature importance'].plot(kind='bar',facecolor='#AA0000')
    ax.patch.set_facecolor('#FFFFFF')
    ax.spines['bottom'].set_color('#CCCCCC')
    ax.spines['bottom'].set_linewidth(1)
    ax.spines['left'].set_color('#CCCCCC')
    ax.spines['left'].set_linewidth(1)
    plt.title('Feature importances', fontsize=24)
# plot ROC with AUC and CI
def plot_ROC(summary_df):

    # plot ROC / PRC for linear model with genes removed
    fig = plt.figure(dpi=400)
    ax = fig.add_subplot(1, 1, 1)
    # df = pd.read_csv('../data/models_analysis/SVM_linear_chip_SS_summary.csv')
    for ind, row in summary_df.iterrows():
        fpr = to_array(row['fpr'])
        tpr = to_array(row['tpr'])

        auc_ci = to_array(row['roc_auc_ci'])

        ax.plot(fpr, tpr,
                label=(', AUC = %0.2f (95%%CI: %0.2f-%0.2f)' 
                      %(row['roc_auc_median'], auc_ci[0], auc_ci[1]))
                )
    ax.plot([0, 1], [0, 1], color='black', lw=2, linestyle=':')
    ax.set_xlim([0.0, 1.0])# plot ROC with AUC and CI
def plot_ROC(summary_df):

    # plot ROC / PRC for linear model with genes removed
    fig = plt.figure(dpi=400)
    ax = fig.add_subplot(1, 1, 1)
    # df = pd.read_csv('../data/models_analysis/SVM_linear_chip_SS_summary.csv')
    for ind, row in summary_df.iterrows():
        fpr = to_array(row['fpr'])
        tpr = to_array(row['tpr'])

        auc_ci = to_array(row['roc_auc_ci'])

        ax.plot(fpr, tpr,
                label=(', AUC = %0.2f (95%%CI: %0.2f-%0.2f)' 
                      %(row['roc_auc_median'], auc_ci[0], auc_ci[1]))
                )
    ax.plot([0, 1], [0, 1], color='black', lw=2, linestyle=':')
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.05])
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    plt.title('ROC')
    plt.legend(loc="lower right")
    # plt.savefig('../Figures/SVM_RBF_allDrugs_genesRemoved_ROC.pdf')
    ax.set_ylim([0.0, 1.05])
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    plt.title('ROC')
    plt.legend(loc="lower right")
    # plt.savefig('../Figures/SVM_RBF_allDrugs_genesRemoved_ROC.pdf')

# plot PRC with AP and CI
def plot_PRC(summary_df):

    fig = plt.figure(dpi=400,)
    ax = fig.add_subplot(1, 1, 1)
#         df = pd.read_csv('../data/analysis/SVM_linear_' + drug_name + '_summary.csv')

    for ind, row in summary_df.iterrows():
        precision = to_array(row['precision'])
        recall = to_array(row['recall'])
    
        ap_ci = to_array(row['AP_ci'])
        ax.step(recall, precision, where='post',
                label=', AP = %0.2f (95%%CI: %0.2f-%0.2f)'
                      %(row['AP_median'], ap_ci[0], ap_ci[1])
                )

    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    ax.set_ylim([0.0, 1.05])
    ax.set_xlim([0.0, 1.0])
    plt.title('PRC')
    plt.legend(loc='lower left')
#     plt.savefig('../Figures/SVM_linear_allDrugs_genesRemoved_PRC.pdf')
    