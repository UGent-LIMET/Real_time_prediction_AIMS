# Title: Software framework for the processing and statistical analysis of multivariate MS-data
# Owner: Laboratory of Integrative Metabolomics (LIMET), Ghent university
# Creator: Dr. Marilyn De Graeve, Volter Paukku
# Maintainer: <limet@ugent.be>
# Script: main.py


##########Configuration##########
#Working directory

#options
PATH_USER = 'xxx/Pipeline_metabolomics/' 

CODE_AUTORUN = 'run code in ternminal automatically'
CODE_DEVELOPMENT = 'run code manually in Rstudio for development locally'

## Adjustments
PATH = PATH_USER

CODE_RUN_MODE = CODE_AUTORUN

#
#####################


##########Sources##########
#source experiment name
import sys
path_Py_scripts = PATH + 'Py_scripts'


if CODE_RUN_MODE == CODE_AUTORUN:
    ## Recognize projectname from python3 commando
    try:
    # argum0 is name of script, 2nd arg=1 is the experiment name
        EXPERIMENT = sys.argv[1]
        name_project = EXPERIMENT
        
        #Source functions
        sys.path.append(path_Py_scripts)
        import data_loading
        import data_statistics
        import data_converting
        import data_writing
        import data_plot
        #import must be in same folder as main.py, write name script wo .py and wo ""
        #does only work when call from terminal, call code always with script.funct
        
        print(data_loading.test_import(name_project)) #test
        
    # test if there is one argument: if not, return an error
    except IndexError:
        raise SystemExit("Your projectname must be supplied")
    
    # test if there is one argument: if not, return an error
    if len(sys.argv) > 2:
        SystemExit("Only one projectname must be supplied")
               
    
if CODE_RUN_MODE == CODE_DEVELOPMENT:
    #run code from configaltenative.py
    #run code from other sources, eg data_functionsALL.py
    name_project = EXPERIMENT

#print(name_project) #test


## data_paths
path_data_in = PATH + 'Data/input/' + name_project + '/'
path_data_out = PATH + 'Data/output/' + name_project + '/' #do not create folder, must be present already
path_data_in_bio = path_data_in + 'bio/'
import pathlib
pathlib.Path(path_data_out).mkdir(parents=True, exist_ok=True)
path_data_out_scans = path_data_out + 'scans/'

## databases paths
path_databases = PATH + 'Databases/'
#not enough space within py box, leave here


## load variables from configalternative
#https://stackoverflow.com/questions/44491442/how-to-import-another-python-script-py-into-main-python-file#44491698
import sys
sys.path.append(path_data_in)
from configalternative import *

print(test) #test correct load config


#
#####################


##########Py Pipeline - Main##########
from datetime import datetime
startTime = datetime.now()
print("Py pipeline - start!")

#Project name:
#!!! todo use assert (name_project == EXPERIMENT) + raise Exception('your experiment name....')
if name_project != EXPERIMENT:
    print("your EXPERIMENT name is incorrect")
    quit()

print(name_project)


#Run selected modules
if RUN_PART_REIMS_TRAINING == RUN_CODE:
    sys.path.append(path_Py_scripts)
    from REIMS_fingerprint_training import * #libraries, funcs
    import REIMS_fingerprint_training
    #run body code module in func
    REIMS_fingerprint_training.run_part_REIMS_training(EXPERIMENT, 
                                                       path_data_in, 
                                                       path_data_out, 
                                                       path_data_out_scans,
                                                       SM,
                                                       TARGET,
                                                       SOURCE_X,
                                                       M3,
                                                       M2,
                                                       VM,
                                                       VARIABLEMETADATA_EXTERN,
                                                       COLLUMN_NR_START_SAMPLES,
                                                       ML_METHOD,
                                                       SVM_GRID,
                                                       RF_GRID,
                                                       NB,
                                                       RF,
                                                       RF_REGRESSION,
                                                       ML_MODEL,
                                                       DONT_SAVE_MODEL,
                                                       SAVE_MODEL
                                                       ) 

    
if RUN_PART_REIMS_PREDICTION == RUN_CODE:
    sys.path.append(path_Py_scripts)
    from REIMS_fingerprint_prediction import * #libraries, funcs
    import REIMS_fingerprint_prediction
    #run body code module in func
    REIMS_fingerprint_prediction.run_part_REIMS_prediction(EXPERIMENT, 
                                                       path_data_in, 
                                                       path_data_out, 
                                                       ML_MODELNAME,
                                                       SOURCE_X,
                                                       M3,
                                                       M2,
                                                       VM,
                                                       VARIABLEMETADATA_EXTERN,
                                                       COLLUMN_NR_START_SAMPLES
                                                       ) 


if RUN_PART_REIMS_EVALUATION == RUN_CODE:
    sys.path.append(path_Py_scripts)
    from REIMS_fingerprint_evaluation import * #libraries, funcs
    import REIMS_fingerprint_evaluation
    #run body code module in func
    REIMS_fingerprint_evaluation.run_part_REIMS_evaluation(EXPERIMENT, 
                                                       path_data_in, 
                                                       path_data_out, 
                                                       SM,
                                                       TARGET,
                                                       Y_PRED
                                                       ) 
 

print(name_project)
print("Py pipeline - done!")
print(datetime.now() - startTime)
#
####################

