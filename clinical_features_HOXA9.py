import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from numpy import mean, var, median
from scipy import stats
from collections import Counter


def find_gene_index(gene_list,gene):
    j = [i for i,x in enumerate(gene_list) if x == gene]
    return j

def find_patients_index(patients, p):
    j = [i for i,x in enumerate(patients) if x == p]
    return j[0]

filename = "log_modified_LAML_TPM.csv"
filename2 = "patients.txt"
filename3 = "FAB.txt"
filename4 = "sex.txt"
filename5 = "age.txt"
filename6 = "BM_blasts.txt"
filename7 = "WBC.txt"
#filename = "modified_raw_counts.csv"
data = pd.read_csv(filename)
patients = pd.read_csv(filename2)
FAB = pd.read_csv(filename3)
sex = pd.read_csv(filename4)
age = pd.read_csv(filename5)
blasts = pd.read_csv(filename6)
WBC = pd.read_csv(filename7)

gene = 'HOXA9'
gene_list = data['Hybridization REF']


# find the index of HOXA9 in the data
i_HOXA9 = find_gene_index(gene_list, gene)
HOXA9_exp = data.iloc[i_HOXA9,2:]

peak1_indexes = [i+2 for i,x in enumerate(HOXA9_exp.values[0]) if x <= 1 and x >= 0.005] # +1 due to the first gene columns we removed +1 due to index shift
peak2_indexes = [i+2 for i,x in enumerate(HOXA9_exp.values[0]) if x <= 5.5 and x >= 4]
peak1_patients = data.iloc[:,peak1_indexes] 
peak2_patients = data.iloc[:,peak2_indexes] 
 
print(len(peak1_indexes))
print(len(peak2_indexes))
 
 
peak1_patients.index = gene_list
peak2_patients.index = gene_list


# only keep the patient number
low_patients_names = [item.split('-')[2] for item in peak1_patients.columns]
high_patients_names = [item.split('-')[2] for item in peak2_patients.columns]


index_cohort1 =[find_patients_index(patients['patients'],int(item)) for item in low_patients_names]
index_cohort2 =[find_patients_index(patients['patients'],int(item)) for item in high_patients_names]
sex_list_cohort1 = sex['sex'][index_cohort1].values.tolist()
sex_list_cohort2 = sex['sex'][index_cohort2].values.tolist()

age_list_cohort1 = age['age'][index_cohort1].values.tolist()
age_list_cohort2 = age['age'][index_cohort2].values.tolist()
FAB_list_cohort1 = FAB['FAB'][index_cohort1].values.tolist()
FAB_list_cohort2 = FAB['FAB'][index_cohort2].values.tolist()




#age_class_low = ['old' if x > median(age_list_cohort1) else 'young' for x in age_list_cohort1]
#age_class_high = ['old' if x > median(age_list_cohort2) else 'young' for x in age_list_cohort2]

print("Age\n")
print(stats.shapiro(age_list_cohort1)) # normal distribution (p>0.05)
print(stats.shapiro(age_list_cohort2)) # NOT
print(stats.levene(age_list_cohort1,age_list_cohort2)) # SAME VARIANCE
print(stats.ttest_ind(age_list_cohort1,age_list_cohort2)) # CANT APPLY T TEST
# multiply pvalue by 2 to have two sided test
print(stats.mannwhitneyu(age_list_cohort1,age_list_cohort2)) # significant

WBC_list_cohort1 = WBC['WBC'][index_cohort1].values.tolist()
WBC_list_cohort2 = WBC['WBC'][index_cohort2].values.tolist()

print("WBC\n")
print(stats.shapiro(WBC_list_cohort1)) # NOT
print(stats.shapiro(WBC_list_cohort2)) # NOT
print(stats.levene(WBC_list_cohort1,WBC_list_cohort2)) # NOT SAME VARIANCE
print(stats.ttest_ind(WBC_list_cohort1,WBC_list_cohort2)) # CANT APPLY T TEST
print(stats.mannwhitneyu(WBC_list_cohort1,WBC_list_cohort2)) # significant

blasts_list_cohort1 = blasts['blasts'][index_cohort1].values.tolist()
blasts_list_cohort2 = blasts['blasts'][index_cohort2].values.tolist()

print("blasts")
print(stats.shapiro(blasts_list_cohort1)) # maybe
print(stats.shapiro(blasts_list_cohort2)) # NOT
print(stats.levene(blasts_list_cohort1,blasts_list_cohort2)) # SAME VARIANCE
print(stats.ttest_ind(blasts_list_cohort1,blasts_list_cohort2)) # CANT APPLY T TEST
print(stats.mannwhitneyu(blasts_list_cohort1,blasts_list_cohort2)) # pvalue x 2 = 0.05 - significant?


d_low = {'age': age_list_cohort1, 'FAB': FAB_list_cohort1, 'blast': blasts_list_cohort1, 'WBC': WBC_list_cohort1, 'HOXA9': ['Low']*len(age_list_cohort1)}
low = pd.DataFrame(data=d_low)
d_high = {'age': age_list_cohort2, 'FAB': FAB_list_cohort2, 'blast': blasts_list_cohort2, 'WBC': WBC_list_cohort2, 'HOXA9': ['High']*len(age_list_cohort2)}
high = pd.DataFrame(data=d_high)

age_FAB =  pd.concat([low,high])
age_FAB = age_FAB.sort_values(by=['FAB'])



Hoxa9_sorted = age_FAB.sort_values(by=['HOXA9'],ascending=False)
plt.subplot(131)
ax1 = sns.boxplot(x='HOXA9',y='age',data=Hoxa9_sorted)
ax1.set_xticklabels(ax1.get_xticklabels(),rotation=30)
ax1.set_xlabel('Age\n ($p=0.0009$)',fontsize=17)
ax1.xaxis.set_label_position('top') 
ax1.set_ylabel('')

plt.subplot(132)
ax2 = sns.boxplot(x='HOXA9',y='WBC',data=Hoxa9_sorted)
ax2.set_xticklabels(ax2.get_xticklabels(),rotation=30)
ax2.set_xlabel('WBC\n ($p=0.0003$)',fontsize=17)
ax2.xaxis.set_label_position('top') 
ax2.set_ylabel('')

plt.subplot(133)
ax3 = sns.boxplot(x='HOXA9',y='blast',data=Hoxa9_sorted)
ax3.set_xticklabels(ax3.get_xticklabels(),rotation=30)
ax3.set_xlabel('%Blast (BM)\n ($p=0.054$)',fontsize=17)
ax3.xaxis.set_label_position('top') 
ax3.set_ylabel('')
plt.show()
