# Title: Software framework for the processing and statistical analysis of multivariate MS-data
# Owner: Laboratory of Integrative Metabolomics (LIMET), Ghent university
# Creator: Dr. Marilyn De Graeve
# Maintainer: <limet@ugent.be>
# Script: functions 


## data_statistics functions

def normalize_with_tic(test):
    #TIC or total ion count normalisation accord to literature
    #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3124646/
    #this normalization will be based on the sum of all intensity values in the mass spectrum (i.e., TIC)
    #= intensity compID 1 / sum(intenstity alls compIDs) for 1 sample/spectrum
    averageTIC_before = test.mean().mean()
    
    #devide by sum per sample
    normalized = test.transpose().sum() #iterate over rows =1 (samples), but saves results as column!!! 
    test = test.apply(lambda num : num / normalized.transpose())   #needs to be transposed for correct tic-norm
      
    #must be re-scaled to same order of scale (afrondingsfouten bij kleine getallen tijdens berek PCA, daarom all getallen maal zelfde factor)
    averageTIC_after = test.mean().mean()
    factor = averageTIC_before/averageTIC_after
    test = test.apply(lambda num : num * factor)

    return(test)


def log_transformation(test): 
    # Log transform: (1+x, base=exp(1))
    #natural log of value+1 so works accurately for negative vanues -1<x<0
    test = test.apply(lambda num : np.log(num + 1))
 
    return(test)


def Pareto_scaling(test):
    # Scaling: pareto: (observation - column mean)/sqrt(stdev column)
    #almost synonym standard scaler (and synonym autoscaler) maar met de wortel SD ipv SD
    test = test.apply(lambda num : num - num.mean())
    test = test.apply(lambda num : num / num.std() ** 0.5)
      
    return(test)


def prediction_Xarray(loaded_model, X_unknown):    
    #list samplenames for output
    samplenames = list(X_unknown.index)
    prediction_output = pd.DataFrame(samplenames, columns=['SampleName'])
    
    ## if predict 1 sample, via array    
    if len(X_unknown) == 1:
        #change to array because 1 sample = 1D array with single feature
        #Reshape your data either using array.reshape(-1, 1) if your data has a single feature or array.reshape(1, -1) if it contains a single sample.
        X_unknown = X_unknown.values.reshape(1,-1)   

    ## predition of unknown X array with predictive ML model
    X_prediction = loaded_model.predict(X_unknown)
    #print("predicted classification number of samples " + str(samplenames) + " are: " + str(X_prediction)) #show only last
    
    ##pred_proba also add
    #X_proba = loaded_model.predict_proba(X_unknown)
    #X_proba = pd.DataFrame(X_proba) #colname from 0, 1, ... 
    
    #save in df
    prediction_output['Predicted_classification'] = X_prediction
    #prediction_output = pd.concat([prediction_output.reset_index(drop=True), X_proba], axis=1)

    return prediction_output



