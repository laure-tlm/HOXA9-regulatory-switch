import pandas as pd
import os
import numpy as np

# filtering and processing of the TPM data

def find_patient(data,patient): # give index of the patient id column -> will be necessary to find TET2 mutated patients
    new_patient_id = [id.split('-')[2] for id in data.columns[1:]]
    if patient not in new_patient_id:
        print('Patient not found\n')
        return
    else:
        i = new_patient_id.index(patient) + 2 # +1 due to absence of first column
    return i

def find_patient_modified(data,patient): # give index of the patient id column -> will be necessary to find TET2 mutated patients
    new_patient_id = [id.split('-')[2] for id in data.columns[2:]]
    if patient not in new_patient_id:
        print('Patient not found\n')
        return
    else:
        i = new_patient_id.index(patient) + 1 # +1 due to absence of first column
    return i

def remove_unknown_genes(data): # remove all the rows with unknown genes and modify the gene id
    new_gene_list = [gene.split('|',1)[0] for gene in data['Hybridization REF']]
    new_gene_list[16301] = 'SLC35E2bis' # to avoid replicates    
    data['Hybridization REF'] = new_gene_list
    indexes = [i for i,x in enumerate(new_gene_list) if x == '?']
    data = data.drop(indexes)
    return data

def low_exp_genes(data): # give index of genes with at least 80 gene exp lower than 1 TPM
    low_exp = [] #index of genes that have a weird low expression
    for i in range(1,len(data)-1):
        low = sum(x<1 for x in data.iloc[i][1:])
        if low>50:
            low_exp.append(i)
    return low_exp

def log_data(data):
    l_data = data.loc[:, data.columns != 'Hybridization REF'].astype(np.float64).applymap(lambda x: np.log(x+1))
    l_data.insert(loc=0,column='Hybridization REF',value=data['Hybridization REF'])
    return l_data

 # open files 
filename = "LAML_TPM.csv"
TPM_data = pd.read_csv(filename)
 
 # remove unknown genes
TPM_data = remove_unknown_genes(TPM_data)
TPM_data.index = range(len(TPM_data))
 
 # rewrite the patient id (maybe not necessary)
new_patients_id = [id.split('-01T',1)[0] for id in TPM_data.columns[1:]]
TPM_data.columns = ['Hybridization REF']+new_patients_id
 
 # remove genes with low gene expression 
low_exp_genes_list = low_exp_genes(TPM_data)
TPM_data = TPM_data.drop(low_exp_genes_list) # remove low expressed genes from our data
TPM_data.index = range(len(TPM_data))
log_TPM_data = log_data(TPM_data)
 
 
out_name = "log_modified_LAML_TPM.csv"
log_TPM_data.to_csv(out_name)


