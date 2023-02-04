
#%% Libraries
import copy 
import itertools
import numpy as np


#%% Datafly algorithm


def datafly(database,functions,k): 

    ''' database : a p*n array where n is the number of quasi-identifier and p is the number of patients
        functions : list of n function for generalization
        k : k-anonymization constraint
        
        output : an anonymised table according to the Datafly algorithm. Some subjects may have been deleted but the structure of database is retained  ''' 
    
    n=len(functions)
    p=len(database)
    
    # function returning the nb of value a QI takes
    def nb_unique_elements(i,array):
        return len(np.unique(array[:,i]))
        
    # function returning how many different patients do not satisfy k-anonymity in the given array
    # if indice is true, it returns the indices of the patients who do no satisfy k-anonymity
    def nb_k_diff_patient(array,indice=False):
        res=0
        indices = []
        for j in range(len(array)):
            if array.count(array[j]) < k:
                res+=1
                if indice:
                    indices.append(j)
        if indice:
            return indices
        return res       
            
    # initialization of the number of different patients in the table
    s_diff = nb_k_diff_patient(database)
    res = copy.deepcopy(database)
    
    while s_diff > k:
        
        # get the QI which contains the most different values
        max_qi = nb_unique_elements(0,res)
        i_max = 0
        for i in range(1,n):
            if nb_unique_elements(i,res) > max_qi:
                i_max=i
                max_qi=nb_unique_elements(i,res)
        
        # The values of this QI are generalized
        for j in range(p):
            res[j][i_max] = res[i_max](res[j][i_max],1) 
            # The 1 indicates that the IQ is generalized to the higher level 
        
        s_diff = nb_k_diff_patient(res)
            
    # On supprime les sujets qui ne respectent pas le k-anonymat
    I = nb_k_diff_patient(res, indice=True)
    for i in I:
        res[i] = ["*" for j in range(n)]
    
    return res