"""
# Title: Software framework for the processing and statistical analysis of multivariate MS-data
# Owner: Laboratory of Integrative Metabolomics (LIMET), Ghent university
# Creator: Dr. Marilyn De Graeve
# Maintainer: <limet@ugent.be>
# Script: Part III cont.: REIMS evaluation prediction 

"""


##########R Pipeline - Part 3: REIMS fingerprint evaluation##########


## load libraries
import time
import pandas as pd
import numpy as np
#https://scikit-learn.org/stable/modules/classes.html#module-sklearn.metrics
import sklearn.metrics as skm
from sklearn.metrics import confusion_matrix
from sklearn.metrics import ConfusionMatrixDisplay
import matplotlib
import matplotlib.pyplot as plt

#Source functions
import data_loading
import data_statistics
import data_converting
import data_writing
import data_plot
#var = 'goeval'
#print(data_loading.test_import(var)) #test

  

#if import in main: need all variables into funct for this to work
#####start func#####
def run_part_REIMS_evaluation(
    EXPERIMENT,
    path_data_in,
    path_data_out,
    SM,
    TARGET,
    Y_PRED,
):

    
    #print(time.time())
    start_time = time.time()
    print("R pipeline - Part III cont.: REIMS fingerprint evaluation - start!")
    # Part III: REIMS fingerprint evaluation
    
      
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
        
     
    
    #load predicted labels from output  
    pred_handle = path_data_out + Y_PRED
    prediction_output = pd.read_csv(pred_handle, sep='\t')
    
    #make y_pred
    y_pred = prediction_output[['SampleName', 'Predicted_classification']]
    y_pred = pd.DataFrame(y_pred)
    y_pred = y_pred.rename(columns={"SampleName": "Index", "Predicted_classification": "Predicted_classification"})
    y_pred = y_pred.set_index('Index') 
    #y_pred.head(2)
    
    
    
    ## only keep non-unique indices
    s1 = pd.merge(y_pred, y, left_index=True, right_index=True) #index is samplename
    #s1.head()
    print("(y_pred.shape), (y.shape) before merge into (s1.shape) are:")
    print(y_pred.shape, y.shape, s1.shape)
    
   
    #pull x,y from each other again
    y_pred = s1.loc[:, s1.columns[0]] #only retain column = "Predicted_classification"
    y = s1.loc[:, s1.columns[-1]] #only retain "Classification"
    print("(y_pred.shape), (y.shape) after merge into (s1.shape) are:")
    print(y_pred.shape, y.shape, s1.shape)
    
    
    
    ### Evaluation ML model (in depth)
    y_test = y
    
    
    ##for manual easy comparison pred + actual
    filename_summary_Y_labels = path_data_out + EXPERIMENT + '_' + ML_MODELNAME[:-4] +'_evaluation_table.txt' 
    df_colnames = s1.columns
    s1.to_csv(filename_summary_Y_labels, header=df_colnames, index=True, sep='\t', mode='w')
    
    ## param write
    #y_pred = classifier.predict(X_test)
    #print(skm.classification_report(y_test,y_pred))
    classification_report = skm.classification_report(y_test,y_pred, output_dict=True)
    classification_report = pd.DataFrame(classification_report).transpose()   
    filename_prediction_output = path_data_out + 'classification_report_' + TARGET + '.txt'
    classification_report.to_csv(filename_prediction_output, header=True, index=True, sep='\t', mode='w')
    
    
    ## confusion matrix print + write + image + plot
    print(confusion_matrix(y_test,y_pred)) 
    
    ##cm write nice
    #add fake data to make correct cm misclassification table (with all members of (m)comp using labels, square matrix)
    y_test_labels = y_test.unique()
    y_pred_labels = y_pred.unique()
    labels = np.append(y_test_labels, y_pred_labels)
    labels = np.unique(labels)
    len(labels)
    
    d = {'Predicted_classication': labels, 'Classification': labels}
    fake_df = pd.DataFrame(d)
    fake_df  
    y_test = pd.concat([y_test, fake_df["Classification"]])
    y_pred = pd.concat([y_pred, fake_df["Predicted_classication"]])
    
    #make cm matrix square
    cm_w_labels = pd.crosstab(y_test,y_pred)   #, margins=True, margins_name = 'Total', dropna=False)   #https://datagy.io/pandas-crosstab/
    
    #rm fake data after achieving square matrix shape
    correct_values = np.diagonal(cm_w_labels) - 1 #1 added for each label in (m)comp
    np.fill_diagonal(cm_w_labels.values, correct_values) 

    #eval
    total = cm_w_labels.sum().sum()
    print("total passes:")
    print(np.diagonal(cm_w_labels).sum())
    print("total fails:")
    print(total - np.diagonal(cm_w_labels).sum())
    print("accuracy:")
    acc = round(np.diagonal(cm_w_labels).sum() / total *100,2)
    print(acc)
    
    #members col
    cm_w_labels.loc[:,'Members'] = cm_w_labels.sum(axis=1) #com members=total
    
    #total row
    cm_w_labels.loc['Total',:] = cm_w_labels.sum(axis=0) #row totalcm_w_labels = pd.crosstab(y_test, y_pred) #row true, col pred

    #correct% col
    i = 0
    cm_w_labels.loc[:,'Correct (%)'] = "" #prep, col=7
    col_nr_correct = cm_w_labels.columns.get_loc('Correct (%)')
    for elem in np.diagonal(cm_w_labels):
        #print(elem)
        #print(i)
        if cm_w_labels.iloc[i,(col_nr_correct-1)] == 0:   #if zero members, otherwise #div/0!
            cm_w_labels.iloc[i,col_nr_correct] = "Not calculated" 
        else:
            cm_w_labels.iloc[i,col_nr_correct] = round(cm_w_labels.iloc[i,i] / cm_w_labels.iloc[i,(col_nr_correct-1)]*100 ,2)  #col correct% = col diag / members col from same row
        i += 1
    
    #accuracy add in last cell
    cm_w_labels.iloc[-1 ,-1] = acc 

    #write
    filename_cm = path_data_out + EXPERIMENT + '_' + ML_MODELNAME[:-4] +'_misclassification_table.txt' 
    cm_w_labels.to_csv(filename_cm, header=True, index=True, sep='\t', mode='w')
    

    
    #image for interactive view (if incomplete: incorr labels)
    cm = confusion_matrix(y_test,y_pred)
    disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels = labels)
    disp.plot(cmap=plt.cm.Blues) 
  
    

    #save cm plot
    #https://matplotlib.org/stable/gallery/images_contours_and_fields/image_annotated_heatmap.html
    plt.figure()
    #plt.title('Confusion Matrix')
    plt.imshow(cm, cmap=plt.cm.Blues)
    threshold_color_text = (cm.max()+cm.min())/2
    for i in range(len(labels)):
       for j in range(len(labels)):  
           if cm[i, j] <= threshold_color_text:
               plt.text(j, i, cm[i, j],
                              ha="center", va="center", color="darkblue")
           else:
               plt.text(j, i, cm[i, j],
                              ha="center", va="center", color="white")
    plt.colorbar()
    plt.xticks(ticks=np.arange(len(labels)), labels=labels)
    plt.yticks(ticks=np.arange(len(labels)), labels=labels)
    plt.xlabel('Predicted label')
    plt.ylabel('True label')
    #plt.tight_layout()
    plt.savefig(path_data_out + 'cm_' + TARGET + '.png', dpi=300)
    
       
    '''
    import matplotlib.pyplot as plt
    from sklearn.metrics import plot_confusion_matrix
    
    #if model still loaded
    classifier = loaded_model
    
    plot_confusion_matrix(classifier, X_test, y_test, cmap=plt.cm.Blues, normalize=None)
    plt.savefig(path_data_out + 'cm_' + TARGET + '.png', dpi=300)
    
    plot_confusion_matrix(classifier, X_test, y_test, cmap=plt.cm.Blues, normalize='true')
    plt.savefig(path_data_out + 'cm_percentage_' + TARGET + '.png', dpi=300)
    '''
    
    
    
    print("R pipeline - Part III cont.: REIMS fingerprint evaluation - done!")
    #print(time.time())
    end_time = time.time()
    print(str(round((end_time - start_time) / 60.0, 2)) + " min")
    #
    #####################
#
#####end func#####     

