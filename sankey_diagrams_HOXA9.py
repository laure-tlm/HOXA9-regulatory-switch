import plotly
import plotly.plotly as py

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from numpy import mean, var
from scipy import stats
from collections import Counter


def find_gene_index(gene_list,gene):
    j = [i for i,x in enumerate(gene_list) if x == gene]
    return j

def find_patients_index(patients, p):
    j = [i for i,x in enumerate(patients) if x == p]
    return j[0]

def how_many_goes(c1,c2):
    j = [x for x in c1 if x in c2]
    return (len(j))

filename = "log_modified_LAML_TPM.csv"
data = pd.read_csv(filename)

filename2 = "patient_description.csv"
patient_description = pd.read_csv(filename2)

patient_description = patient_description[:-1]


gene = 'HOXA9'


gene_list = data['Hybridization REF']

data.index = gene_list
data = data.drop(columns='Hybridization REF').transpose()

data2 = data[[gene]].sort_values(by=[gene])


peak1 = data2[  (data2[gene] >= 0.005)  & (data2[gene] <= 1)]
peak2 = data2[  (data2[gene] >= 3)  & (data2[gene] <= 6)]


peak1_patients_names = [item.split('-')[2] for item in peak1.index]
peak2_patients_names = [item.split('-')[2] for item in peak2.index]

peak1_FAB = [ patient_description['FAB'][i] for i,x in enumerate(patient_description['TCGA Patient ID']) if str(int(x)) in peak1_patients_names]
peak2_FAB = [ patient_description['FAB'][i] for i,x in enumerate(patient_description['TCGA Patient ID']) if str(int(x)) in peak2_patients_names]

peak1_classification = [patient_description['Molecular Classification'][i] for i,x in enumerate(patient_description['TCGA Patient ID']) if str(int(x)) in peak1_patients_names]
peak2_classification = [patient_description['Molecular Classification'][i] for i,x in enumerate(patient_description['TCGA Patient ID']) if str(int(x)) in peak2_patients_names]

peak1_d = {'Molecular classification': peak1_classification,  'FAB': peak1_FAB,gene: ['low']*len(peak1_FAB)}
peak1 = pd.DataFrame(data=peak1_d)
peak2_d = {'Molecular classification': peak2_classification,  'FAB': peak2_FAB,gene: ['High']*len(peak2_FAB)}
peak2 = pd.DataFrame(data=peak2_d)

#print(Counter(peak1['Molecular classification']))
print(Counter(peak2['Molecular classification']))

#################### PLOT SANKEY DIAGRAM #########

Molecular = dict(
    type='sankey',
    node = dict(
      pad = 5,
      thickness = 15,
      line = dict(
        color = "black",
        width = 2
      ),
      label = ["Low peak HOXA9","High peak HOXA9","Normal Karyotype","Intermediate Risk Cytogenic Abn.","PML-RARA","CBF translocations","Complex Cytogenetics","MLL/NUP98 translocations","Poor Risk Cytogenic Abn."],
      color = ["blue", "orange", "grey","grey","blue","blue","orange","orange","orange"]
    ),
    link = dict(
      source = [0,0,0,0,1,1,1,1,1],
      target = [2,3,4,5,2,3,6,7,8],
      value = [8,2,11,9,62,10,21,11,7]
  ))

layout =  dict(
    #title = "Sankey Diagram for HOXA9/APP and HOXA9/APP/IGSF10 FAB clusters",
    font = dict(
      size = 20
    )
)

FAB = dict(
    type='sankey',
    node = dict(
      pad = 5,
      thickness = 15,
      line = dict(
        color = "black",
        width = 2
      ),
      label = ["Low peak HOXA9","High peak HOXA9","M0","M1","M2","M3","M4","M5"],
      color = ["blue", "orange", "orange", "grey", "grey","blue","grey","orange"]
    ),
    link = dict(
      source = [0,0,0,0,1,1,1,1,1],
      target = [3,4,5,6,2,3,4,6,7],
      value = [2,12,11,6,15,34,20,22,18]
  ))

layout =  dict(
    #title = "Sankey Diagram for HOXA9/APP and HOXA9/APP/IGSF10 FAB clusters",
    font = dict(
      size = 20
    )
)


fig = dict(data=[Molecular], layout=layout)
#py.iplot(fig, validate=False)
plotly.offline.plot(fig, validate=False)
