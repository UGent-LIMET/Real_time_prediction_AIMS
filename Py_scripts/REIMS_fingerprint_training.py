"""
# Title: Software framework for the processing and statistical analysis of multivariate MS-data
# Owner: Laboratory of Integrative Metabolomics (LIMET), Ghent university
# Creator: Dr. Marilyn De Graeve
# Maintainer: <limet@ugent.be>
# Script: Part III: REIMS fingerprint training

"""


##########R Pipeline - Part 3: REIMS fingerprint training ##########


## load libraries
import time
import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn import svm
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import GaussianNB
from sklearn.pipeline import Pipeline, FeatureUnion
from sklego.meta import EstimatorTransformer
import sklearn.model_selection

import random

import numpy as np
import matplotlib.pyplot as plt
#https://scikit-learn.org/stable/modules/classes.html#module-sklearn.metrics
import sklearn.metrics as skm
from sklearn.metrics import confusion_matrix
from sklearn.metrics import plot_confusion_matrix

import pickle
import os



#Source functions
import data_loading
import data_statistics
import data_converting
import data_writing
import data_plot
#var = 'gotrain'
#print(data_loading.test_import(var)) #test



#if import in main: need all variables into funct for this to work
#####start func#####
def run_part_REIMS_training(
    EXPERIMENT,
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
    SAVE_MODEL,
):

    #print(time.time())
    start_time = time.time()
    print("R pipeline - Part III: training REIMS fingerprints - start!")
    # Part III: REIMS fingerprint training


    ## Yarray
    #load preY
    sampleMetadata_handle = path_data_in + SM
    sampleMetadata = pd.read_csv(sampleMetadata_handle, sep='\t')
    #sampleMetadata.shape
    #sampleMetadata.head()
    
    #keep rows with info groeps
    comp = sampleMetadata.dropna(subset = [TARGET]) #keep 1,2,3... or -1/+1 not nan
    
    #make Y
    y = comp[['SampleName', TARGET]]
    y = pd.DataFrame(y)
    y = y.rename(columns={"SampleName": "Index", TARGET: "Classification"})
    y = y.set_index('Index')
    #y.head(2) 
    #y.index #with '...' ipv Int64Index
    
    
    
    ## perform fitting according to source
    if (SOURCE_X == M3 or SOURCE_X == VM):
        
        ## load depending source
        if SOURCE_X == M3:
            ## load and prepare X array
            X = data_loading.load_Xarray_from_M3()
        if SOURCE_X == VM:
            ## load and prepare X array
            vm_handle = path_data_in + VARIABLEMETADATA_EXTERN
            X = data_loading.load_Xarray_from_VM(vm_handle, COLLUMN_NR_START_SAMPLES)
            #X.index #make sure no "X..." in name #with '...'
            #X.head(2)
    
    
        ##perform TIC normalization by default
        #test = X.iloc[:10,0:10]
        #test = normalize_with_tic(test)
        X = data_statistics.normalize_with_tic(X)
        
  
    
        ## make train/test
        ## only keep non-unique indices
        s1 = pd.merge(X, y, left_index=True, right_index=True) #index is samplename
        #s1.head()
        print("(X.shape), (y.shape) before merge into (s1.shape) are:")
        print(X.shape, y.shape, s1.shape)
    
       
    
        #pull x,y from each other again
        X = s1.drop(s1.columns[-1],axis=1) #drop last column = "Classification"
        y = s1.loc[:, s1.columns[-1]] #only retain "Classification"
        print("(X.shape), (y.shape) after merge into (s1.shape) are:")
        print(X.shape, y.shape, s1.shape)
    
    
        ##split X and y
        X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=1)
        
        #check sizes dfs
        print("X_train.shape, y_train.shape, X_test.shape, y_test.shape are:")
        print(X_train.shape, y_train.shape, X_test.shape, y_test.shape)
    
    
        ## balance train set
        if BALANCE == DONT_BALANCE:
            print('hi')
            #do nothing
        if BALANCE == BALANCE_SAME_SIZE:
            s1 = pd.merge(X_train, y_train, left_index=True, right_index=True)
            if "Comparison" in TARGET:
                #balance if pairwise comparison
                if "Multi" not in TARGET:
                    s1 = data_converting.BalanceDF(s1) #using -1 and +1
                    print("(s1.shape) from test after balancing is:")
                    print(s1.shape)
                #balance if pairwise comparison
                if "MultipleComparison" in TARGET:
                    s1 = data_converting.BalanceDF_multi(s1)
                    print("(s1.shape) from test after balancing is:")
                    print(s1.shape)
            X_train = s1.drop(s1.columns[-1],axis=1) #drop last column = "Classification"
            y_train = s1.loc[:, s1.columns[-1]] #only retain "Classification"      
            print("X_train.shape, y_train.shape after balancing are:")
            print(X_train.shape, y_train.shape)
            
        if BALANCE == BALANCE_RESAMPLING:
            s1 = pd.merge(X_train, y_train, left_index=True, right_index=True)
            if "Comparison" in TARGET:
                #balance if pairwise comparison
                if "Multi" not in TARGET:
                    s1 = data_converting.BalanceDF_by_resampling(s1)
                    print("(s1.shape) from test after balancing is:")
                    print(s1.shape)
                #balance if pairwise comparison
                if "MultipleComparison" in TARGET:
                    s1 = data_converting.BalanceDF_by_resampling(s1)
                    print("(s1.shape) from test after balancing is:")
                    print(s1.shape)
            X_train = s1.drop(s1.columns[-1],axis=1) #drop last column = "Classification"
            y_train = s1.loc[:, s1.columns[-1]] #only retain "Classification"      
            print("X_train.shape, y_train.shape after balancing are:")
            print(X_train.shape, y_train.shape)

    
        ##perform log transformation and Pareto scaling (not needed, SS done in pipe instead)
        #test = X.iloc[:10,0:10]
        #test = log_transformation(test)
        #test = Pareto_scaling(test)
        
        #X_train = log_transformation(X_train)
        #X_test = log_transformation(X_test)
        #X_train = Pareto_scaling(X_train)
        #X_test = Pareto_scaling(X_test)
        
            
        ##train ML model according to chosen ML_METHOD
    
        ################################
        if ML_METHOD == SVM_GRID:
            steps = [('ss', StandardScaler()), ('SVM', svm.SVC())]
            pipeline = Pipeline(steps) # define the pipeline object.
            
            parameteres = {'SVM__kernel': ['linear', 'poly', 'rbf', 'sigmoid', 'precomputed'], 
                           'SVM__C':[0.001,0.1,10,100,10e5], 
                           'SVM__gamma':[0.1,0.01]}
            #=> if: Invalid parameter kernel for estimator Pipeline(steps=[('ss', StandardScaler()), ('SVM', SVC())]). 
            #Check the list of available parameters with `estimator.get_params().keys()
            #pipeline.get_params().keys()
            
            #gridsearch
            #grid = sklearn.model_selection.RandomizedSearchCV(pipeline, param_distributions=parameteres, cv=5, scoring="accuracy")
            grid = sklearn.model_selection.GridSearchCV(pipeline, param_grid=parameteres, cv=5, scoring="accuracy")
            
            #run the gridsearch
            grid.fit(X_train, y_train)
            
            ## same name use for model
            classifier = grid 
            #pipeline #name before export model
            
            print("best model parameters were:")
            print(classifier.best_params_)
        #################################
    
        
        ################################
        if ML_METHOD == RF_GRID:
            steps = [('ss', StandardScaler()), ('RF', RandomForestClassifier())]
            pipeline = Pipeline(steps) # define the pipeline object.
             
            parameteres = {'RF__n_estimators':[10, 25, 50, 100, 250, 500], 
                           'RF__max_depth':[15, 30, 75, 125, None], 
                           'RF__min_samples_leaf': [1, 3], 
                           'RF__max_leaf_nodes': [10, 20, 35, 50, None]}
            #=> if: Invalid parameter kernel for estimator Pipeline(steps=[('ss', StandardScaler()), ('SVM', SVC())]). 
            #Check the list of available parameters with `estimator.get_params().keys()
            #pipeline.get_params().keys()
            
            #gridsearch
            grid = sklearn.model_selection.RandomizedSearchCV(pipeline, param_distributions=parameteres, cv=5, scoring="accuracy")
            
            #run the gridsearch
            grid.fit(X_train, y_train)
            
            ## same name use for model
            classifier = grid 
            #pipeline #name before export model
            
            print("best model parameters were:")
            print(classifier.best_params_)
        #################################

    
        ################################
        if ML_METHOD == LDA_GRID: 
           
            steps = [('ss', StandardScaler()), ('LDA', LinearDiscriminantAnalysis())]
            pipeline = Pipeline(steps) # define the pipeline object.
             
            parameteres = {'LDA__n_components':[5, 10, 15, 20, 25, 30]}
            #todo: literature, which other most import to test?
            #pipeline.get_params().keys()
            
            #gridsearch
            grid = sklearn.model_selection.GridSearchCV(pipeline, param_grid=parameteres, cv=5, scoring="accuracy")
            
            #run the gridsearch
            grid.fit(X_train, y_train)
            
            ## same name use for model
            classifier = grid 
            #pipeline #name before export model
            
            print("best model parameters were:")
            print(classifier.best_params_)
        #################################


        ################################
        if ML_METHOD == NB:
            steps = [('ss', StandardScaler()), ('nb', GaussianNB())]
            pipeline = Pipeline(steps) # define the pipeline object.
             
            parameteres = {}
            #=> if: Invalid parameter kernel for estimator Pipeline(steps=[('ss', StandardScaler()), ('SVM', SVC())]). 
            #Check the list of available parameters with `estimator.get_params().keys()
            #pipeline.get_params().keys()
            
            #gridsearch
            grid = sklearn.model_selection.GridSearchCV(pipeline, param_grid=parameteres, cv=5)
            
            #run the gridsearch
            grid.fit(X_train, y_train)
               
            ## same name use for model
            classifier = grid 
            #pipeline #name before export model
            
            print("best model parameters were:")
            print(classifier.best_params_)
        #################################
        
        
        #################################
        if ML_METHOD == RF_REGRESSION:
            #try wo grid => no, no relaible, min 5x CV
            #regr = RandomForestRegressor() #nothing fixed, score with default max_depth=2, random_state=0
            #regr.fit(X_train, y_train)   
            #classifier = regr
            #print("model parameters were:")
            #print(classifier.get_params)
            
            steps = [('RF', RandomForestRegressor())]
            pipeline = Pipeline(steps) # define the pipeline object.
             
            parameteres = {}
            #=> if: Invalid parameter kernel for estimator Pipeline(steps=[('ss', StandardScaler()), ('SVM', SVC())]). 
            #Check the list of available parameters with `estimator.get_params().keys()
            #pipeline.get_params().keys()
            
            #gridsearch
            grid = sklearn.model_selection.RandomizedSearchCV(pipeline, param_distributions=parameteres, cv=5)
            
            #run the gridsearch
            grid.fit(X_train, y_train)
            
            ## same name use for model
            classifier = grid 
            #pipeline #name before export model
            
            print("best model parameters were:")
            print(classifier.best_params_)
        #################################
             
        
        ## evaluate model 
        print("accuracy of model:")
        print(classifier.score(X_test,y_test)) 
            # Default metric is R2 for regression, which can be accessed by score()
            # Default metric is accuracy for clasifci
        print("for target:")
        print(TARGET)
        print('for ML method:')
        print(ML_METHOD)
        
    
    
    if SOURCE_X == M2:
        ## load and prepare X array + train models (grid, all) per sample => all in one loop
        
        print('STOP: M2 option deprived in code. Use M3 or VM')
        quit()
        
        #list all M2 files in scans folder, no subset for now
        filenames_X = []
        directory_list = os.listdir(path_data_out_scans)
        
        #random order list with filenames
        random.seed(4)
        random.shuffle(directory_list)
        
        for filename in directory_list:
            #print (filename) #all files, folders
            if "matrix2_" in filename:
                filenames_X.append(path_data_out_scans + filename)
            
        
        #amount of files to devide per 100 
        #nr_filenames_X = len(filenames_X) 
        #nr_groups_to_train = int(round(nr_filenames_X /100, 0))
        nr_groups_to_train = 1
        groups_to_train = np.array_split(filenames_X, nr_groups_to_train)
        #print(groups_to_train[0])
        
        
        #iter per subset of samples: in this loop: split, merge, train model via partial fit
        for group in groups_to_train:
            #group = groups_to_train[0]
            #print(group)
            
            ## make train/test
            X_train = pd.DataFrame()
            X_test = pd.DataFrame()
            filenames_X_train = []
            filenames_X_test = []
            i = 0
            for filename in group:
                #print(filename) #all files, folders
                #filename = filenames_X[5]
                 
                #split test-train per group
                if i % 3 == 0: 
                    #1/3th of data is test set, rest in train
                    #print(i)
                    filenames_X_test.append(filename)
                    X_test_1s = data_loading.load_Xarray_from_M2_filename(filename)
                    X_test = X_test.append(X_test_1s)    
                else:
                    filenames_X_train.append(filename)
                    X_train_1s = data_loading.load_Xarray_from_M2_filename(filename)
                    X_train = X_train.append(X_train_1s)
                i = i + 1
    
                
            
            ## only keep non-unique indices (train)
            s1 = pd.merge(X_train, y, left_index=True, right_index=True) #index is samplename
            #s1.head()
            print("(X_train.shape), (y.shape) before merge into (s1.shape) are:")
            print(X_train.shape, y.shape, s1.shape)
            
            #balance if pairwise comparison
            if "Multi" not in TARGET:
                s1 = data_converting.BalanceDF(s1)
                print("(s1.shape) from train after balancing is:")
                print(s1.shape)
            
            #pull x,y from each other again
            X_train = s1.drop(s1.columns[-1],axis=1) #drop last column = "Classification"
            y_train = s1.loc[:, s1.columns[-1]] #only retain "Classification"
            print("(X_train.shape), (y_train.shape) after merge into (s1.shape) are:")
            print(X_train.shape, y_train.shape, s1.shape)
               
            
            ## only keep non-unique indices (test)
            s1 = pd.merge(X_test, y, left_index=True, right_index=True) #index is samplename
            #s1.head()
            print("(X_test.shape), (y.shape) before merge into (s1.shape) are:")
            print(X_test.shape, y.shape, s1.shape)
            
            #balance if pairwise comparison
            if "Multi" not in TARGET:
                s1 = data_converting.BalanceDF(s1)
                print("(s1.shape) from test after balancing is:")
                print(s1.shape)
                
            #pull x,y from each other again
            X_test = s1.drop(s1.columns[-1],axis=1) #drop last column = "Classification"
            y_test = s1.loc[:, s1.columns[-1]] #only retain "Classification"
            print("(X_test.shape), (y_test.shape) after merge into (s1.shape) are:")
            print(X_test.shape, y_test.shape, s1.shape)
        
        
            #check sizes dfs
            print("X_train.shape, y_train.shape, X_test.shape, y_test.shape are:")
            print(X_train.shape, y_train.shape, X_test.shape, y_test.shape)
               
            
            ################################
            steps = [('ss', StandardScaler()), ('RF', RandomForestClassifier())]
            pipeline = Pipeline(steps) # define the pipeline object.
             
            parameteres = {'RF__n_estimators':[10, 25, 50, 100, 250, 500], 
                           'RF__max_depth':[15, 30, 75, 125, None], 
                           'RF__min_samples_leaf': [1, 3], 
                           'RF__max_leaf_nodes': [10, 20, 35, 50, None]}
            #=> if: Invalid parameter kernel for estimator Pipeline(steps=[('ss', StandardScaler()), ('SVM', SVC())]). 
            #Check the list of available parameters with `estimator.get_params().keys()
            #pipeline.get_params().keys()
            
            #gridsearch
            grid = sklearn.model_selection.RandomizedSearchCV(pipeline, param_distributions=parameteres, cv=5)
            
            #run the gridsearch
            grid.fit(X_train, y_train)
            #################################
            
            
        ### after iter all from all M2 (all groups per loop)     
        ## same name use for model
        classifier = grid 
        #pipeline #name before export model
        
        
        ## evaluate model 
        print("accuracy of model and param are:")
        print(classifier.score(X_test,y_test))
        print(classifier.best_params_)
    
    
    
    ### Evaluation ML model
    if ML_METHOD != RF_REGRESSION:
        ### Evaluation ML classification model (in depth)
        #works only for classification
        ## param write
        y_pred = classifier.predict(X_test)
        #print(skm.classification_report(y_test,y_pred))
        classification_report = skm.classification_report(y_test,y_pred, output_dict=True)
        classification_report = pd.DataFrame(classification_report).transpose()   
        filename_prediction_output = path_data_out + 'classification_report_' + TARGET + '.txt'
        classification_report.to_csv(filename_prediction_output, header=True, index=True, sep='\t', mode='w')
        
            
        ## confusion matrix print + save pngs
        print(confusion_matrix(y_test,y_pred))
        
        plot_confusion_matrix(classifier, X_test, y_test, cmap=plt.cm.Blues, normalize=None)
        plt.savefig(path_data_out + 'cm_' + TARGET + '.png', dpi=300)
        
        plot_confusion_matrix(classifier, X_test, y_test, cmap=plt.cm.Blues, normalize='true')
        plt.savefig(path_data_out + 'cm_percentage_' + TARGET + '.png', dpi=300)
    
    
    if ML_METHOD == RF_REGRESSION:
        ### Evaluation ML regression model (in depth)
        R2_score = classifier.score(X_test,y_test)
        
         ## ! all in one file, call funct train 100x but in 1doc. must be present in output (write wo header info)
        # Open a file with access mode 'a'
        file_object = open(path_data_out + EXPERIMENT + '_Predicted_regressions.txt' , 'a')
        # Append 'hello' at the end of file
        #file_object.write("\n")
        file_object.write(TARGET)
        file_object.write('\t')
        file_object.write(str(R2_score))
        file_object.write('\n')
        # Close the file
        file_object.close()
        

    
    if ML_MODEL == SAVE_MODEL:
        # save the model to disk
        timestr = time.strftime("%H%M%S")
        filename = path_data_out + EXPERIMENT + '_MLmodel_' + TARGET + '_' + timestr + '.pkl' #adjust the good ones
        pickle.dump(classifier, open(filename, 'wb'))   #adjust the good ones



    print("R pipeline - Part III: REIMS fingerprint training - done!")
    #print(time.time())
    end_time = time.time()
    print(str(round((end_time - start_time) / 60.0, 2)) + " min")
    
#
#####end func#####



#
#####################

