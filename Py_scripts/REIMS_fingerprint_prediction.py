"""
# Title: Software framework for the processing and statistical analysis of multivariate MS-data
# Owner: Laboratory of Integrative Metabolomics (LIMET), Ghent university
# Creator: Dr. Marilyn De Graeve
# Maintainer: <limet@ugent.be>
# Script: Part III cont.: REIMS fingerprint prediction 

"""


##########R Pipeline - Part 3: REIMS fingerprint prediction##########


## load libraries
import time
import pickle
import pandas as pd
import os
    
    
    
#Source functions
import data_loading
import data_statistics
import data_converting
import data_writing
import data_plot
#var = 'gopred'
#print(data_loading.test_import(var)) #test



#if import in main: need all variables into funct for this to work
#####start func#####
def run_part_REIMS_prediction(
    EXPERIMENT,
    path_data_in,
    path_data_out,
    ML_MODELNAME,
    SOURCE_X,
    M3,
    M2,
    VM,
    VARIABLEMETADATA_EXTERN,
    COLLUMN_NR_START_SAMPLES,
):
    

    #print(time.time())
    start_time = time.time()
    print("R pipeline - Part III cont.: REIMS fingerprint prediction - start!")
    # Part III: REIMS fingerprint prediction
    
    
    ## load the saved ML model from disk
    filename = path_data_in + ML_MODELNAME
    loaded_model = pickle.load(open(filename, 'rb'))
    #result = loaded_model.score(X_test1, y_test1)
    #print(result)
    
    
    
    ## perform prediction according to source
    if SOURCE_X == M3:
        ## load and prepare X array
        X_unknown = data_loading.load_Xarray_from_M3()
      
        ## prediction unknown Xarray
        prediction_output = data_statistics.prediction_Xarray(loaded_model, X_unknown)
        
        #write output to file
        filename_predicted_Y_labels = path_data_out + EXPERIMENT + '_' + ML_MODELNAME[:-4] + '_Predicted_classifications.txt' 
        df_colnames = prediction_output.columns
        prediction_output.to_csv(filename_predicted_Y_labels, header=df_colnames, index=False, sep='\t', mode='w')
    
    
    if SOURCE_X == VM:
        ## load and prepare X array
        vm_handle = path_data_in + VARIABLEMETADATA_EXTERN
        X_unknown = data_loading.load_Xarray_from_VM(vm_handle, COLLUMN_NR_START_SAMPLES)
        
        ##perform TIC normalization by default
        X_unknown = data_converting.normalize_with_tic(X_unknown)
        
        ## prediction unknown Xarray
        prediction_output = data_statistics.prediction_Xarray(loaded_model, X_unknown)
        
        #write output to file
        filename_predicted_Y_labels = path_data_out + EXPERIMENT + '_' + ML_MODELNAME[:-4] +'_Predicted_classifications.txt' 
        df_colnames = prediction_output.columns
        prediction_output.to_csv(filename_predicted_Y_labels, header=df_colnames, index=False, sep='\t', mode='w')
    
    
    if SOURCE_X == M2:
        ## load and prepare X array + predict per sample
        #append results in prediction_output
        prediction_output = data_loading.load_predict_Xarray_from_M2() 
        
        #write output to file
        filename_predicted_Y_labels = path_data_out + EXPERIMENT + '_' + ML_MODELNAME[:-4] + '_Predicted_classifications.txt' 
        df_colnames = prediction_output.columns
        file_exists = os.path.isfile(filename_predicted_Y_labels)
        if file_exists:
            prediction_output.to_csv(filename_predicted_Y_labels, header=False, index=False, sep='\t', mode='a')
        else:
            prediction_output.to_csv(filename_predicted_Y_labels, header=df_colnames, index=False, sep='\t', mode='a')
    
        
    
    
    print("R pipeline - Part III cont.: REIMS fingerprint prediction - done!")
    #print(time.time())
    end_time = time.time()
    print(str(round((end_time - start_time) / 60.0, 2)) + " min")
    #
    #####################
#
#####end func#####    

