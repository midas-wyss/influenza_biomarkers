# -*- coding: utf-8 -*-
"""
Created on Sun Jun 16 15:34:22 2019

@author: malcantar
"""

import pandas as pd
import numpy as np

from sklearn.metrics import roc_curve, auc
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from sklearn.model_selection import train_test_split

from eli5.sklearn import PermutationImportance
from utils import one_hot_encode
from utils import summary_stats

from tqdm import tqdm
import pickle
from utils import get_weights
import mygene
from biothings_client import get_client





# does SVM predictions and outputs important metrics as a dataframe
#
def pred_SVM(X_train, X_test, y_train, y_test, model):
    model.fit(X_train, y_train)
    try: #for SVM
        y_score = model.predict(X_test)
        y_score = model.decision_function(X_test)
    except AttributeError: #for RF
        y_score = model.predict_proba(X_test)[:, 1]
    fpr, tpr, _ = roc_curve(y_test, y_score)
    roc_auc = auc(fpr, tpr)
    
    averagePrecision = average_precision_score(y_test, y_score)
    precision, recall, _ = precision_recall_curve(y_test, y_score)
    
    modelAccuracy = model.score(X_test, y_test)
    results_df = pd.DataFrame({'fpr': [fpr], 'tpr': [tpr], 
                               'roc_auc': roc_auc,
                               'precision': [precision], 
                               'recall': [recall],
                               'AP': averagePrecision,
                               'accuracy': modelAccuracy})
    perm = PermutationImportance(model).fit(X_test, y_test)
    return results_df, perm

# function for grabbing top features for linear SVM
def feat_ML(X_train, y_train, model):
    
    model.fit(X_train, y_train)
    try: #for SVM 
        weights = pd.DataFrame([model.named_steps["clf"].coef_[0]], 
                                columns=X_train.columns[model.named_steps["fs"].get_support()])
    except AttributeError: #for RF
        try:
            weights = pd.DataFrame([model.feature_importances_],
                                   columns=X_train.columns)
        except AttributeError:
            weights = None 
    return weights

# run linear SVM after training 
def runLinearSVM(count_data_zScore_df, sample_metadata_df, numSimulations, training_size_percent,
                  MODEL_LOC, MODEL_STATS_LOC,
                 save=False):
    
    simulation_iterations = numSimulations

    sim_weights_all_df = pd.DataFrame()
    all_weights_ens = []
    all_results_ens = []

    virus_control_labels = one_hot_encode(count_data_zScore_df, sample_metadata_df)
    # loading model from pickle file
    with open(MODEL_LOC + '.pkl', 'rb') as model_pickle:
        model = pickle.load(model_pickle)

    for sim in tqdm(range(simulation_iterations)):
        X_train, X_test, y_train, y_test = train_test_split(count_data_zScore_df, virus_control_labels, 
                                                            train_size = training_size_percent, 
                                                            stratify=virus_control_labels, random_state = sim)

        # predict using SVM
        all_results, perm_import = pred_SVM(X_train, X_test, y_train, y_test, model)
        all_results['seed'] = sim
        all_results_ens.append(all_results)

        # feature weights
        sim_weights_df = pd.DataFrame([model.named_steps["clf"].coef_[0]], 
                                columns=count_data_zScore_df.columns[model.named_steps["fs"].get_support()])
        # concatenate feature weights
        sim_weights_all_df = pd.concat([sim_weights_all_df, sim_weights_df], ignore_index=True)

        # grab feature weights
        all_weights = feat_ML(X_train, y_train, model)
        if all_weights is not None:
            all_weights['seed'] = sim
            all_weights_ens.append(all_weights)
    
    # concatenating feature weights from all simulations and calculating summary statistics 
    all_results_ens = pd.concat(all_results_ens, sort=False).set_index('seed')
    all_results_sum = summary_stats(all_results_ens)


    # exporting feature weights and other summary stats from simulations 
    if all_weights_ens:
        all_weights_ens = pd.concat(all_weights_ens, sort=False).set_index('seed')
        #Normalize weights per row to get relative importances
        all_weights_ens = get_weights(all_weights_ens)
        all_weights_ens.to_csv(MODEL_STATS_LOC+'_all_weights.csv')

    if save:
        all_results_sum.to_csv(MODEL_STATS_LOC +'_summary.csv')
        all_results_ens.to_csv(MODEL_STATS_LOC  +'_ensResults.csv')
        sim_weights_all_df.to_csv(MODEL_STATS_LOC +'_weights_EachSim.csv')
#     return(perm_import)
        mg = mygene.MyGeneInfo()
        symbol_to_entrez = pd.read_csv('../data/gene_lists/gene_symbol_to_entrez.csv', index_col = 0)
        symbol_to_entrez_dict = dict(zip(symbol_to_entrez['SYMBOL'], 
                                       symbol_to_entrez['ENTREZID']))
        feature_importances_df = all_weights_ens.copy()
        feature_importances_annotation_df = feature_importances_df.copy()
        for gene in feature_importances_df.index:
            try:
                entrez_ID_temp = symbol_to_entrez_dict[gene]
                argument_temp = 'entrezgene: ' + str(entrez_ID_temp)
                name_temp = mg.query(argument_temp)['hits'][0]['name']
                feature_importances_annotation_df.loc[gene,'annotation'] = name_temp
            except IndexError:
                 feature_importances_annotation_df.loc[gene,'annotation'] = 'None'
        feature_importances_annotation_df.to_csv(MODEL_STATS_LOC  +'_all_weights_annotation.csv')

def evaluate_RF(model, test_features, test_labels):
    predictions = model.predict(test_features)
    accuracy = 1-np.sum(abs(predictions-test_labels)) / len(predictions)
    return accuracy

# run linear SVM after training 
def run_RF(count_data_zScore_df, sample_metadata_df, numSimulations, training_size_percent,
                  MODEL_LOC, MODEL_STATS_LOC,
                 save=False):
    
    simulation_iterations = numSimulations

    sim_weights_all_df = pd.DataFrame()
    all_weights_ens = []
    all_results_ens = []

    virus_control_labels = one_hot_encode(count_data_zScore_df, sample_metadata_df)
    # loading model from pickle file
    with open(MODEL_LOC + '.pkl', 'rb') as model_pickle:
        model = pickle.load(model_pickle)

    for sim in tqdm(range(simulation_iterations)):
        X_train, X_test, y_train, y_test = train_test_split(count_data_zScore_df, virus_control_labels, 
                                                            train_size = training_size_percent, 
                                                            stratify=virus_control_labels, random_state = sim)

        # predict using SVM
        all_results, perm_import = pred_SVM(X_train, X_test, y_train, y_test, model)
        all_results['seed'] = sim
        all_results_ens.append(all_results)

        sim_weights_df = pd.DataFrame(data=model.feature_importances_, 
                        columns = ['Importance'], index=list(count_data_zScore_df.columns))
        sim_weights_df = sim_weights_df.copy().T
        sim_weights_all_df = pd.concat([sim_weights_all_df, sim_weights_df], ignore_index=True)

         # grab feature weights
        all_weights = feat_ML(X_train, y_train, model)
        if all_weights is not None:
            all_weights['seed'] = sim
            all_weights_ens.append(all_weights)
    
    # concatenating feature weights from all simulations and calculating summary statistics 
    all_results_ens = pd.concat(all_results_ens, sort=False).set_index('seed')
    all_results_sum = summary_stats(all_results_ens)

    # exporting feature weights and other summary stats from simulations 
    if all_weights_ens:
        all_weights_ens = pd.concat(all_weights_ens, sort=False).set_index('seed')
        #Normalize weights per row to get relative importances
        all_weights_ens = get_weights(all_weights_ens)
        all_weights_ens.to_csv(MODEL_STATS_LOC+'_all_weights.csv')
    if save:
        all_results_sum.to_csv(MODEL_STATS_LOC +'_summary.csv')
        all_results_ens.to_csv(MODEL_STATS_LOC  +'_ensResults.csv')
        sim_weights_all_df.to_csv(MODEL_STATS_LOC +'_weights_EachSim.csv')
#      return(perm_import)
        mg = mygene.MyGeneInfo()
        symbol_to_entrez = pd.read_csv('../data/gene_lists/gene_symbol_to_entrez.csv', index_col = 0)
        symbol_to_entrez_dict = dict(zip(symbol_to_entrez['SYMBOL'], 
                                       symbol_to_entrez['ENTREZID']))
        feature_importances_df = all_weights_ens.copy()
        feature_importances_annotation_df = feature_importances_df.copy()
        for gene in feature_importances_df.index:
            try:
                entrez_ID_temp = symbol_to_entrez_dict[gene]
                argument_temp = 'entrezgene: ' + str(entrez_ID_temp)
                name_temp = mg.query(argument_temp)['hits'][0]['name']
                feature_importances_annotation_df.loc[gene,'annotation'] = name_temp
            except IndexError:
                 feature_importances_annotation_df.loc[gene,'annotation'] = 'None'
        feature_importances_annotation_df.to_csv(MODEL_STATS_LOC  +'_all_weights_annotation.csv')
