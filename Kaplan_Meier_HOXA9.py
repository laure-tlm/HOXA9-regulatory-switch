import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from numpy import mean, var
from scipy import stats
from matplotlib import rc
from lifelines import KaplanMeierFitter

# python program to plot the OS difference between young/old and high/low HOXA9 patients

def find_gene_index(gene_list,gene):
    j = [i for i,x in enumerate(gene_list) if x == gene]
    return j

def find_patients_index(patients, p):
    j = [i for i,x in enumerate(patients) if x == p]
    return j[0]

filename = "log_modified_LAML_TPM.csv"
filename2 = "laml_tcga_clinical_data.tsv" # from David download - cbioPortal

data = pd.read_csv(filename)
patient_description = pd.read_csv(filename2,sep='\t')

gene_list = data['Hybridization REF']

# find the index of HOXA9 in the data
i_HOXA9 = find_gene_index(gene_list, "HOXA9")
HOXA9_exp = data.iloc[i_HOXA9,2:]

# select patients that have HOXA9 expression in the peaks
peak1_indexes = [i+2 for i,x in enumerate(HOXA9_exp.values[0]) if x <= 1 and x >= 0.005] # +1 due to the first gene columns we removed +1 due to index shift
peak2_indexes = [i+2 for i,x in enumerate(HOXA9_exp.values[0]) if x <= 5.5 and x >= 4]

# 32 patients for low and 80 for high

peak1_patients = data.iloc[:,peak1_indexes].columns
peak2_patients = data.iloc[:,peak2_indexes] .columns

# only keep the patient number
peak1_patients = [item.split('-')[2] for item in peak1_patients]
peak2_patients = [item.split('-')[2] for item in peak2_patients]

patient2 = patient_description['Patient ID']
patient2 = [item.split('-')[2] for item in patient2]

p1_indexes =  [i for i,x in enumerate(patient2) if x in peak1_patients]
p2_indexes =  [i for i,x in enumerate(patient2) if x in peak2_patients]

p1_vital = patient_description["Patient's Vital Status"][p1_indexes]
p2_vital = patient_description["Patient's Vital Status"][p2_indexes]

p1_vital2 = [0 if item == "Alive" else 1 for item in p1_vital]
p2_vital2 = [0 if item == "Alive" else 1 for item in p2_vital]

p1_OS = patient_description["Overall Survival (Months)"][p1_indexes]
p2_OS = patient_description["Overall Survival (Months)"][p2_indexes]

p1_age = patient_description["Diagnosis Age"][p1_indexes]
p2_age = patient_description["Diagnosis Age"][p2_indexes]

p1 = {'OS':p1_OS, 'vital':p1_vital, 'vital2':p1_vital2, 'age':p1_age}
p2 = {'OS':p2_OS, 'vital':p2_vital, 'vital2':p2_vital2, 'age':p2_age}

p1 = pd.DataFrame(data=p1)
p2 = pd.DataFrame(data=p2)

p1["HOXA9"] = "Low"
p2["HOXA9"] = "High"

p1["age cat"] = ["old" if i > 60 else "young" for i in p1["age"]]
p2["age cat"] = ["old" if i > 60 else "young" for i in p2["age"]]

p1["cluster"] = p1["HOXA9"] +"/"+ p1["age cat"]
p2["cluster"] = p2["HOXA9"] +"/"+ p2["age cat"]

cohorts_patients =  pd.concat([p1,p2])
cohorts_patients.to_csv('HOXA9_survival_age.csv')




