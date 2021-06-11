import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from numpy import mean, var
from scipy import stats
from matplotlib import rc
from lifelines import KaplanMeierFitter

# python program to plot the OS difference between M2 HOXA9 low and M2 high HOXA9

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

M2_low_indexes =  [i for i,x in enumerate(patient2) if x in peak1_patients and patient_description['FAB'][i] == 'M2']
M2_high_indexes =  [i for i,x in enumerate(patient2) if x in peak2_patients and patient_description['FAB'][i] == 'M2']

M4_low_indexes =  [i for i,x in enumerate(patient2) if x in peak1_patients and patient_description['FAB'][i] == 'M4']
M4_high_indexes =  [i for i,x in enumerate(patient2) if x in peak2_patients and patient_description['FAB'][i] == 'M4']


M2_low_vital = patient_description["Patient's Vital Status"][M2_low_indexes]
M2_high_vital = patient_description["Patient's Vital Status"][M2_high_indexes ]
M4_low_vital = patient_description["Patient's Vital Status"][M4_low_indexes]
M4_high_vital = patient_description["Patient's Vital Status"][M4_high_indexes ]

M2_low_vital2 = [0 if item == "Alive" else 1 for item in M2_low_vital]
M2_high_vital2 = [0 if item == "Alive" else 1 for item in M2_high_vital]
M4_low_vital2 = [0 if item == "Alive" else 1 for item in M4_low_vital]
M4_high_vital2 = [0 if item == "Alive" else 1 for item in M4_high_vital]

M2_low_OS = patient_description["Overall Survival (Months)"][M2_low_indexes]
M2_high_OS = patient_description["Overall Survival (Months)"][M2_high_indexes]
M4_low_OS = patient_description["Overall Survival (Months)"][M4_low_indexes]
M4_high_OS = patient_description["Overall Survival (Months)"][M4_high_indexes]

M2_low_tab = {'OS':M2_low_OS, 'vital':M2_low_vital, 'vital2':M2_low_vital2}
M2_high_tab = {'OS':M2_high_OS, 'vital':M2_high_vital, 'vital2':M2_high_vital2}
M4_low_tab = {'OS':M4_low_OS, 'vital':M4_low_vital, 'vital2':M4_low_vital2}
M4_high_tab = {'OS':M4_high_OS, 'vital':M4_high_vital, 'vital2':M4_high_vital2}

M2_low_tab = pd.DataFrame(data=M2_low_tab)
M2_high_tab = pd.DataFrame(data=M2_high_tab)
M4_low_tab = pd.DataFrame(data=M4_low_tab)
M4_high_tab = pd.DataFrame(data=M4_high_tab)

M2_low_tab["HOXA9"] = "Low"
M2_high_tab["HOXA9"] = "High"
M4_low_tab["HOXA9"] = "Low"
M4_high_tab["HOXA9"] = "High"

M2_patients =  pd.concat([M2_low_tab,M2_high_tab])
M4_patients =  pd.concat([M4_low_tab,M4_high_tab])

M2_patients.to_csv('M2_HOXA9_survival.csv')
M4_patients.to_csv('M4_HOXA9_survival.csv')



