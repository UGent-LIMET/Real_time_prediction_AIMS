# Title: Software framework for the processing and statistical analysis of multivariate MS-data
# Owner: Laboratory of Integrative Metabolomics (LIMET), Ghent university
# Creator: Dr. Marilyn De Graeve
# Maintainer: <limet@ugent.be>
# Script: Configuration alternative


##########Global_settings##########

## Options
RUN_CODE = 'run this part of the pipeline'
DONT_RUN_CODE = 'skip this part of the pipeline, keeps order: pre-processing - targeted analysis - statistical analysis - annotation'

## Adjustments
#Project name:
EXPERIMENT = 'Real_time_prediction_AIMS' 
USER_COMMENT = "Tutorial comment" #Add info about experiment, eg. explain (Multiple)Comparisons, to include in reports
    
RUN_PART_REIMS_TRAINING = RUN_CODE 
RUN_PART_REIMS_PREDICTION = DONT_RUN_CODE
RUN_PART_REIMS_EVALUATION = DONT_RUN_CODE

#
####################



########## REIMS fingerprint training ##########
if RUN_PART_REIMS_TRAINING == RUN_CODE:
    ## options
    #options SOURCE_X
    M3 = 'takes training information from raw MS-spectra (matrix3_CopmIDs_all), use for small datasets'
    M2 = 'takes training information from raw MS-spectra (individual matrix3_CopmIDs files), use for large datasets'
    VM = 'takes training information from pre-processed feature matrices (variableMetadata), use for prediction whitin the same experiment (samples need to be pre-processed together)'
    
    #balancing dataset
    BALANCE_SAME_SIZE = 'discards the excess of samples to match the smallest group'
    BALANCE_RESAMPLING = 'resample samples to match the largest group'
    DONT_BALANCE = 'keep all samples in model training'

    #options ML_METHOD: algorithm and setting to use
    RF_GRID = 'RF with random search used predifined parameters, optimize using CV'
    SVM_GRID = 'SVM with grid search used predefined parameters, optimize using CV'
    LDA_GRID = 'LDA with grid search used predefined parameters, optimize using CV'
    SVM = 'SVM with default parameters, optimize using CV'
    NB = 'naive bayes'
    RF = 'RF with default parameters'
    RF_REGRESSION = 'RF regression, optimize using CV'

    #options to save trained ML_MODEL for future predictions
    DONT_SAVE_MODEL = 'do not save ML model'
    SAVE_MODEL = 'save ML model'
    
    
    ## Adjustments   
    # file containing target (Y)
    SM = '20200828_REIMS_R_fish_sampleMetadataWO5.txt'
    TARGET = 'MultipleComparison1' #"MultipleComparison1" #'Comparison1' #'1' compid1 column if SM_VM_merged file
    
    # file containing inputdata (X)
    SOURCE_X = VM
    
    #in case VM is chosen as SOURCE_X, give filename
    VARIABLEMETADATA_EXTERN = 'TEST_REIMS_fingerprints_training_variableMetadataWO5.txt' 
    COLLUMN_NR_START_SAMPLES = 20

    #ML algoritm
    ML_METHOD = RF_GRID
    
    #save model for future use (prediction)
    ML_MODEL = SAVE_MODEL

#
####################



########## REIMS fingerprint prediction ##########
if RUN_PART_REIMS_PREDICTION == RUN_CODE:
    ## options
    #options SOURCE_X
    M3 = 'takes training information from raw MS-spectra (matrix3_CopmIDs_all), use for small datasets'
    M2 = 'takes training information from raw MS-spectra (individual matrix3_CopmIDs files), use for large datasets'
    VM = 'takes training information from pre-processed feature matrices (variableMetadata), use for prediction whitin the same experiment (samples need to be pre-processed together)'
    
    
    ## Adjustments    
    ML_MODELNAME = "TEST_REIMS_fingerprints_training_MLmodel_MultipleComparison1_150356.pkl" #must be present in input folder of your experiment (todo ev. 1 fixed path 'DBs' later)
    
    SOURCE_X = VM
    
    #in case VM is chosen as SOURCE_X, give filename
    VARIABLEMETADATA_EXTERN = "TEST_REIMS_fingerprints_prediction_variableMetadataW5.txt"
    COLLUMN_NR_START_SAMPLES = 20

#
####################



########## REIMS fingerprint prediction ##########
if RUN_PART_REIMS_EVALUATION == RUN_CODE:    
    ## Adjustments   
    # file containing target (Y)
    SM = '20200828_REIMS_R_fish_sampleMetadataW5.txt' #"20200828_REIMS_R_fish_sampleMetadata.txt" 
    TARGET = 'MultipleComparison1' #'Comparison1   #"MultipleComparison1" 

    Y_PRED = "TEST_REIMS_fingerprints_prediction_TEST_REIMS_fingerprints_training_MLmodel_MultipleComparison1_150356_Predicted_classications.txt"

#
####################



test = 'goconfig' #test loading config ok

