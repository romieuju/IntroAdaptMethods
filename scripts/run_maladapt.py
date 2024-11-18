#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#    IntroAdapt: joint inference of introgression and selection
#    Copyright (C) 2024  Ghislain Camarata and Jules Romieu, ISEM/CBGP, CNRS/INRAE. 
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.

#Script to classify genetic data with MaLadapt, use feature save in pre-trained model save in "config_file/maladapt_feature_lists/" or change path in feature_file variable.
#The path of the models trained is recorded in the variable maladapt_trained_path (ini file)
#Name of the models trained is recorded in the variable maladapt_features (ini file)
#Scaler for each trained model are save on maladapt_trained_path 
#This script should be run using the maladapt environment


import pandas as pd
from inspect import signature
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import chi2
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve
from sklearn.preprocessing import label_binarize
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import cross_val_score
import pickle
from introadapt import *

#Read project option file 
options = get_project_options(proj_options_file = sys.argv[1])

#Save results folder path in a variable (results_dir/analyse/project/project_sim)
project_dir = options["results_dir"]+"/"+options["analysis"]+"/"+options["project"]

maladapt_trained_path = options["maladapt_trained_path"]
maladapt_features     = options["maladapt_features"]

#load data. For summary stats, column order might be relevant
complete_table = pd.read_csv(project_dir+"/maladapt_tab_"+options["project"]+".csv", sep=",")

#replaces unusable areas (the same was done for training)
complete_table = complete_table.replace(float("inf"), 0) 
complete_table = complete_table.fillna(0)
complete_table = complete_table.replace(str("Nan"), 0)

out_table     = complete_table.copy()#creates a hard copy of the table that will receive the output
all_col_names = list(complete_table.columns.values)

for features in maladapt_features:#different models were trained using two different sets of statistics. stat list and model files should be named following the same structure
    model_path   = maladapt_trained_path+"/"+features+"_model.sav"
    scaler_path  = maladapt_trained_path+"/scaler_"+features+".pkl"
    scaler       = pickle.load(open(scaler_path, "rb"))
    ML           = pickle.load(open(model_path, "rb"))
    feature_file = "config_file/maladapt_feature_lists/"+features+".txt"
    with open(feature_file,"rt") as file:
        current_features = file.readline().split()
    
    unpredictive               = [i for i in all_col_names if i not in current_features]    #creates a list of unused features
    cut_table                  = complete_table.drop(unpredictive, axis=1) #drops all columns not found in the stat list file
    test_X                     = cut_table.values#turns the dataframe into an array
    test_X[:,:]                = scaler.transform(test_X[:,:])#standardize data
    y_score                    = ML.predict_proba(test_X[:,:].astype(float))#applies the model to our data
    out_table["mal_"+features] = y_score[:,1] #adds the score given by this model to the output table
    

out_table.to_csv(project_dir+"/maladapt_out_"+options["project"]+".csv", mode='w', header=True, index=False)
    



