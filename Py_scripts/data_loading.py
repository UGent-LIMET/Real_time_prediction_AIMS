# Title: Software framework for the processing and statistical analysis of multivariate MS-data
# Owner: Laboratory of Integrative Metabolomics (LIMET), Ghent university
# Creator: Dr. Marilyn De Graeve
# Maintainer: <limet@ugent.be>
# Script: functions 


## test main.py
def test_import(name_project):
    return name_project + name_project


## data_loading dependencies
import os
import pandas as pd
import numpy as np

  
## data_loading functions

def load_Xarray_from_M3():
    #load preX
    matrix3_handle = path_data_out_scans + "matrix3_CopmIDs_all.txt"
    matrix3_CopmIDs_all = pd.read_csv(matrix3_handle, sep='\t')
    matrix3_CopmIDs_all = matrix3_CopmIDs_all.set_index("MZ") #mz no column but the index
    
    #rm "X" in beginning samplename(s)
    list_samples = list(matrix3_CopmIDs_all.columns)
    #matrix3_CopmIDs_all.columns[0][1:]
    list_samples = [e[1:] for e in list_samples] #wo "X" at start
    matrix3_CopmIDs_all.columns = list_samples
    
    #matrix3_CopmIDs_all.shape
    #test = matrix3_CopmIDs_all.iloc[0:30, 0:30] #mz as index
    #test
    
    #make X
    X = matrix3_CopmIDs_all.T
    #X.shape
    #test = X.iloc[0:30, 0:30] #samplename as index
    #test
    
    #QUB data diff mz range solve (so set to 0); interpol correct 50-1200 so in this range points present to train
    X = X.fillna(0)

    return X


def load_Xarray_from_VM(vm_handle, COLLUMN_NR_START_SAMPLES):
    matrix3_CopmIDs_all = pd.read_csv(vm_handle, sep='\t')
    matrix3_CopmIDs_all = matrix3_CopmIDs_all.set_index("MZ") #mz no column but the index

    #from metadata to matrix
    #matrix3_CopmIDs_all.columns
    matrix3_CopmIDs_all = matrix3_CopmIDs_all.iloc[:,(COLLUMN_NR_START_SAMPLES-2):] #MZ col eruit in funct load as index
    
    #make X
    X = matrix3_CopmIDs_all.T
    
    #QUB data diff mz range so solve temp
    X = X.fillna(0)

    return X


def load_Xarray_from_M2_filename(filename):
        #load Xarrays for each sample in group
        X = pd.read_csv(filename, sep='\t')
        X = X.set_index("MZ") #mz no column but the index
        
        matrix3_CopmIDs_all = X
        #rm "X" in beginning samplename(s)
        list_samples = list(matrix3_CopmIDs_all.columns)
        list_samples = [e[1:] for e in list_samples] #wo "X" at start
        matrix3_CopmIDs_all.columns = list_samples
        
        #make X
        X = matrix3_CopmIDs_all.T
        
        #QUB data diff mz range solve (so set to 0); interpol correct 50-1200 so in this range points present to train
        X = X.fillna(0)
        
        return(X)


def load_predict_Xarray_from_M2(): 
    #list all M2 files in scans folder, no subset for now
    filenames_X = []
    directory_list = os.listdir(path_data_out_scans)
    for filename in directory_list:
        #print (filename) #all files, folders
        if "matrix2_" in filename:
            filenames_X.append(path_data_out_scans + filename)
        
       
    #load files from subset filenames_X list
    prediction_output = pd.DataFrame() #empty df
    for filename in filenames_X:
        #print(filename)
        #filename = sublist_X[1]
        X = pd.read_csv(filename, sep='\t')
        X = X.set_index("MZ") #mz no column but the index
        
        matrix3_CopmIDs_all = X
        #rm "X" in beginning samplename(s)
        list_samples = list(matrix3_CopmIDs_all.columns)
        list_samples = [e[1:] for e in list_samples] #wo "X" at start
        matrix3_CopmIDs_all.columns = list_samples
        
        #make X
        X = matrix3_CopmIDs_all.T
        
        #QUB data diff mz range solve (so set to 0); interpol correct 50-1200 so in this range points present to train
        X = X.fillna(0)
        
        ## prediction unknown Xarray
        prediction_x = prediction_Xarray(loaded_model, X)
        
        #append from all samples to df
        prediction_output = prediction_output.append(prediction_x)

    return prediction_output
        
