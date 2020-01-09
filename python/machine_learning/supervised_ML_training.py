# -*- coding: utf-8 -*-
"""
Created on Sun Jun 16 15:29:50 2019

@author: malcantar
"""
import pprint as pp
import pandas as pd
from sklearn.linear_model import SGDClassifier
import numpy as np
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import SelectFromModel
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import RepeatedStratifiedKFold
# from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RandomizedSearchCV


# package for saving ML models to disk
import pickle
from utils import one_hot_encode
from supervised_ML_testing import evaluate_RF

# train a linear SVM given a feature matrix and sample metadata
# inputs: normalized count matrix, sample metadata, where you want model to save, where you want training stats
# to save, and whether you want to save anything
# outputs: trained model, training results
def train_linear_SVM(count_data_zScore_df,sample_metadata_df, MODEL_LOC, MODEL_STATS_LOC, save = False):
    
    # one-hot encoding control_virus labels
    optimal_feature_dict = {}
    
    virus_control_labels = one_hot_encode(count_data_zScore_df,sample_metadata_df)

    # splitting data
    X_train, _, y_train, _ = train_test_split(count_data_zScore_df, virus_control_labels, 
                                              train_size = 0.8, stratify=virus_control_labels, random_state=42)

    # cross-validation sets
    cv = RepeatedStratifiedKFold(n_splits=int(len(X_train)/2),n_repeats=2)

    # estimator --> hinge loss makes this an SVM
    clf = SGDClassifier(loss="hinge", penalty="l1", shuffle=True, 
                        class_weight="balanced", alpha=0.0001, max_iter=1000, tol=1e-3) 

    # feature selection
    feat_selection = SelectFromModel(clf, threshold= 0.0)

    # SGD-SVM feature selection
    model = Pipeline([
              ('fs', feat_selection),  # SelectFromModel using SGD-SVM
              ('clf', clf), 
            ])

    # Set the parameter ranges for optimizing grid search cross-validation 
    params = {
        'fs__threshold': np.linspace(0,2, num=40), 
        'clf__alpha': np.logspace(-2, -1, num=3)
            }

    # conducting grid search using CV
    grid = GridSearchCV(model , params, cv=cv)
    grid.fit(X_train, y_train)

    print("The best parameters are %s with a score of %0.2f"
          % (grid.best_params_, grid.best_score_))

    # saving best param to dictionary
    optimal_feature_dict.update({"fs__threshold":list(grid.best_params_.values())[1],
                                               'clf__alpha': list(grid.best_params_.values())[0],
                                               "best_score": grid.best_score_})
    # saving best model in pickle format
    pickle_file = MODEL_LOC + '.pkl'

    if save:
        with open(pickle_file, 'wb') as f:
             pickle.dump(grid.best_estimator_, f)

        # saving results from best model
        CVresults = pd.DataFrame(data=grid.cv_results_, 
                                columns=['param_fs__threshold', 
                                        'param_clf__alpha', 
                                        'mean_test_score'])
        CVresults.to_csv(MODEL_STATS_LOC + '.csv')
# trains RandomForest estimator; see SVM training function for description of inputs and outputs
def train_random_forest_estimator(count_data_zScore_df, sample_metadata_df, 
                                  MODEL_LOC, MODEL_STATS_LOC, save = False):

    virus_control_labels = one_hot_encode(count_data_zScore_df, sample_metadata_df)
    
    # number of trees
    n_estimators = [int(x) for x in np.linspace(start = 10, stop = 1000, num = 20)]
    # features to consider at splits
    max_features = ['auto', 'sqrt']
    # max levels in tree
    max_depth = [int(x) for x in np.linspace(1, 20, num = 10)]
    max_depth.append(None)
    # min number of samples required to split
    min_samples_split = [2, 3, 4, 5]
    # min number of samples required at each leaf node
    min_samples_leaf = [1, 2, 4]
    # boostrap
    bootstrap = [True, False]
    
    # Create the random grid
    random_grid = {'n_estimators': n_estimators,
                   'max_features': max_features,
                   'max_depth': max_depth,
                   'min_samples_split': min_samples_split,
                   'min_samples_leaf': min_samples_leaf,
                   'bootstrap': bootstrap}
    pp.pprint(random_grid)
    # Create a based model
    rf = RandomForestClassifier()
    # Instantiate the grid search model
    grid_search = GridSearchCV(estimator = rf, param_grid = random_grid, 
                              cv = 3, n_jobs = -1, verbose = 2)
    X_train, X_test, y_train, y_test = train_test_split(count_data_zScore_df, virus_control_labels, 
                                                                train_size = 0.8, 
                                                                stratify=virus_control_labels, random_state = 42)
    rf_random = RandomizedSearchCV(estimator = rf, param_distributions = random_grid, n_iter = 100, cv = 3, verbose=2, random_state=42, n_jobs = -1)
    # Fit the random search model
    rf_random.fit(X_train, y_train)
    
    # saving best model in pickle format
    pickle_file = MODEL_LOC + '.pkl'

    best_random = rf_random.best_estimator_
    random_accuracy = evaluate_RF(best_random, X_test, y_test)
    
    if save:
        with open(pickle_file, 'wb') as f:
            pickle.dump(best_random, f)  

    CVresults = pd.DataFrame(data=rf_random.cv_results_)
    CVresults.to_csv(MODEL_STATS_LOC + '.csv')
    print('The best predictor obtain an accuracy of:', str(random_accuracy))
