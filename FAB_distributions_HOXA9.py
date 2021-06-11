import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from numpy import mean, var
from scipy import stats
from matplotlib import rc

# python program to plot the difference of subtype percentage between low and high cohorts of HOXA9 expression in AML

def find_gene_index(gene_list,gene):
    j = [i for i,x in enumerate(gene_list) if x == gene]
    return j

def find_patients_index(patients, p):
    j = [i for i,x in enumerate(patients) if x == p]
    return j[0]

filename = "log_modified_LAML_TPM.csv"
filename2 = "patients.txt"
filename3 = "FAB.txt"
#filename = "modified_raw_counts.csv"
data = pd.read_csv(filename)
patients = pd.read_csv(filename2)
FAB = pd.read_csv(filename3)

gene_list = data['Hybridization REF']

# find the index of APP in the data
i_HOXA9 = find_gene_index(gene_list, "HOXA9")
HOXA9_exp = data.iloc[i_HOXA9,2:]

# select patients that have HOXA9 expression in the peaks
peak1_indexes = [i+2 for i,x in enumerate(HOXA9_exp.values[0]) if x <= 1 and x >= 0.005] # +1 due to the first gene columns we removed +1 due to index shift
peak2_indexes = [i+2 for i,x in enumerate(HOXA9_exp.values[0]) if x <= 5.5 and x >= 4]


# 31 patients for low and 80 for high

peak1_patients = data.iloc[:,peak1_indexes].columns
peak2_patients = data.iloc[:,peak2_indexes] .columns

# only keep the patient number
peak1_patients = [item.split('-')[2] for item in peak1_patients]
peak2_patients = [item.split('-')[2] for item in peak2_patients]

# gives the index of the patients and then its associated FAB
FAB_index_low =[find_patients_index(patients['patients'],int(item)) for item in peak1_patients ]
FAB_list_low = FAB['FAB'][FAB_index_low].values.tolist()
FAB_list_low = [i for i in FAB_list_low if i != 'M7' and i != 'M6' and i != 'nc']

FAB_index_high =[find_patients_index(patients['patients'],int(item)) for item in peak2_patients ]
FAB_list_high = FAB['FAB'][FAB_index_high].values.tolist()
FAB_list_high = [i for i in FAB_list_high if i != 'M7' and i != 'M6' and i != 'nc']

# in order to plot a stacked bar plot
M0_means = [100*FAB_list_low.count('M0')/len(FAB_list_low),100*FAB_list_high.count('M0')/len(FAB_list_high)]
M1_means = [100*FAB_list_low.count('M1')/len(FAB_list_low),100*FAB_list_high.count('M1')/len(FAB_list_high)]
M2_means = [100*FAB_list_low.count('M2')/len(FAB_list_low),100*FAB_list_high.count('M2')/len(FAB_list_high)]
M3_means = [100*FAB_list_low.count('M3')/len(FAB_list_low),100*FAB_list_high.count('M3')/len(FAB_list_high)]
M4_means = [100*FAB_list_low.count('M4')/len(FAB_list_low),100*FAB_list_high.count('M4')/len(FAB_list_high)]
M5_means = [100*FAB_list_low.count('M5')/len(FAB_list_low),100*FAB_list_high.count('M5')/len(FAB_list_high)]

# plot
r = [0,1]
barWidth = 0.3
names = ('Low','High')
# Create M0 bars
p1 = plt.bar(r, M0_means, color='b', edgecolor='white', width=barWidth)
# Create M1 bars
p2 = plt.bar(r, M1_means, bottom = M0_means, color='g', edgecolor='white', width=barWidth)
p3 = plt.bar(r, M2_means, bottom=[i+j for i,j in zip(M0_means,M1_means)], color='r', edgecolor='white', width=barWidth)
p4 = plt.bar(r, M3_means, bottom=[i+j+k for i,j,k in zip(M0_means,M1_means,M2_means)], color='c', edgecolor='white', width=barWidth)
p5 = plt.bar(r, M4_means, bottom=[i+j+k+l for i,j,k,l in zip(M0_means,M1_means,M2_means,M3_means)], color='m', edgecolor='white', width=barWidth)
p6 = plt.bar(r, M5_means, bottom=[i+j+k+l+m for i,j,k,l,m in zip(M0_means,M1_means,M2_means,M3_means,M4_means)], color='y', edgecolor='white', width=barWidth)
 
plt.legend((p1[0],p2[0],p3[0],p4[0],p5[0],p6[0]), ('M0', 'M1','M2', 'M3','M4','M5'))

# Custom x axis
plt.xticks(r, names)
plt.xlabel("HOXA9 expression")
plt.ylabel("FAB percentage")
 
# Show graphic
plt.savefig('FABpc_HOXA9_thres1.png')




