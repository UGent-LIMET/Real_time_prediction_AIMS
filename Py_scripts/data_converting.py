# Title: Software framework for the processing and statistical analysis of multivariate MS-data
# Owner: Laboratory of Integrative Metabolomics (LIMET), Ghent university
# Creator: Dr. Marilyn De Graeve
# Maintainer: <limet@ugent.be>
# Script: functions 


## data_converting functions

def BalanceDF(s1):
    #balance same amount F/NF label (with less 0 than 1)
    zeros_subset = s1.loc[s1["Classification"] == -1, :] 
    number_of_0s = len(zeros_subset)
    #print(number_of_0s)
    
    ones_subset = s1.loc[s1["Classification"] == 1, :]
    number_of_1s = len(ones_subset)
    #print(number_of_1s)
    
    if number_of_1s >= number_of_0s:
        sampled_ones = ones_subset.sample(n=number_of_0s, random_state=1)  
        balanced_s1 = pd.concat([zeros_subset, sampled_ones], ignore_index=True)
    if number_of_0s >= number_of_1s:
        sampled_zeros = zeros_subset.sample(n=number_of_1s, random_state=1)  
        balanced_s1 = pd.concat([ones_subset, sampled_zeros], ignore_index=True)
    
    return balanced_s1


def BalanceDF_multi(s1):
    #balance same amount labels in Mcomp
    comp_nr = s1["Classification"].unique() 
    smallest = len(s1["Classification"])  #smallest nr of elem in Mcomp groups, start high
    balanced_s1 = pd.DataFrame() #start empty
    for elem in comp_nr:
        #print(elem)        
        elem_subset = s1.loc[s1["Classification"] == elem, :] 
        number_of_elem = len(elem_subset)
        #print(number_of_elem)
          
        if number_of_elem < smallest:
            smallest = number_of_elem #get smallest amount
    
    for elem in comp_nr: #keep smallest from each group
        #print(elem)     
        elem_subset = s1.loc[s1["Classification"] == elem, :] 
        sampled_elem = elem_subset.sample(n=smallest, random_state=1)  
        balanced_s1 = pd.concat([balanced_s1, sampled_elem], ignore_index=True)
        #print(len(balanced_s1["Classification"]))    
        
    return balanced_s1
    

def BalanceDF_by_resampling(s1):
    #balance same amount labels in Mcomp, by adding replicates from smaller groups
    comp_nr = s1["Classification"].unique() 
    largest = len(s1["Classification"])  #smallest nr of elem in Mcomp groups, start high
    balanced_s1 = pd.DataFrame() #start empty
    for elem in comp_nr:
        #print(elem)        
        elem_subset = s1.loc[s1["Classification"] == elem, :] 
        number_of_elem = len(elem_subset)
        #print(number_of_elem)
          
        if number_of_elem > smallest:
            largest = number_of_elem #get smallest amount
    
    for elem in comp_nr: #keep smallest from each group
        #print(elem)     
        elem_subset = s1.loc[s1["Classification"] == elem, :] 
        sampled_elem = elem_subset.sample(n=largest, replace=True, random_state=1)  
        balanced_s1 = pd.concat([balanced_s1, sampled_elem], ignore_index=True)
        #print(len(balanced_s1["Classification"]))    
        
    return balanced_s1

