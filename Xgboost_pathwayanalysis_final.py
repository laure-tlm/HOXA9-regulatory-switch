## Runs xgboost on cancer pathways - runs as: python machine_learning_pathway_analysis.py GENEOFINTEREST (eg. TP53)

## Load a ton of possibly unnecessary libraries... ##
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib 
matplotlib.rcParams['figure.figsize'] = [15, 10]
plt.rcParams.update({"font.size":20})
import numpy as np
import glob
import os.path
from os import path
import sys
import shap
from scipy.stats import zscore


from sklearn import model_selection
from xgboost import XGBClassifier
from xgboost import plot_importance
from sklearn.model_selection import train_test_split
#from sklearn.cross_validation import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.metrics import matthews_corrcoef
## Dependency on mygene removed, but data must contain row names that match both your cancer pathways and your gene query

import warnings
warnings.simplefilter(action='ignore')

import pickle


##Generate a dictionary of pathways from the following folder
pathwaylist = glob.glob("./../../Cancer_pathways/*.txt")

pathways = {}
for item in pathwaylist:
    pathwayname = (item.split("/")[-1:])[0].split(".")[0]
    print("Setting up genes for " + pathwayname + " pathway")
    genelist = []
    genefile = open(item, "r")
    for line in genefile.readlines():
        if line.startswith(">"):
            continue
        else:
            gene = (line.split("\t")[0]).split("|")[0]
            genelist.append(gene)
        pathways.update({pathwayname:genelist})
print("Saving pathways as pathways_pickle.p")
pickle.dump(pathways, open("pathways_pickle.p", "wb"))
# Save pathways as a pickle (though this is now super-fast)


## Read in gene name argument
genename = sys.argv[1]
## Find the file containing all the data (in this case this is all GI cancers, and the OCCAMS dataset)
cancerfile = "./modified_LAML_TPM.csv"

## Transform the data around a bit
RNAseq_data = pd.read_csv(cancerfile, index_col = 'Hybridization REF')
RNAseq_data = RNAseq_data.loc[~(RNAseq_data==0).all(axis=1)]
RNAseq_t = RNAseq_data.transpose()
RNAseq_newindex = RNAseq_t.drop_duplicates() # Dropping duplicates of gene names mostly removed genes named "???" 
RNAseq_newindex = RNAseq_newindex.apply(zscore)

Combined_gene_sorted = RNAseq_newindex.sort_values(by = [genename], ascending=False) # Sort data by gene of interest
gene_high = Combined_gene_sorted.head(30)
gene_low = Combined_gene_sorted.tail(30) # Take top and bottom 100 patients
gene_high["gene_status"] = 1
gene_low["gene_status"] = 0
Combined_gene = pd.concat([gene_high, gene_low]) # Result is an array with all genes, for the top and bottom 100 expressers of gene X

## Set up a loop to do the actual machine learning
pathwaystats = []
for pathway, genelist in pathways.items(): # Iterate through each pathway and gene list
    print("Analysing " + pathway)
    print(genelist)
    if not os.path.exists(pathway):
        os.makedirs(pathway) # Make a folder for each pathway
    os.chdir(pathway)
    for i in range(0,100): # Build 100 different models (with a different random seed)
        importance_array = "Importance_%i" % (i+1)
        accuracy_array = "Accuracy_%i" % (i+1)
        shap_array = "SHAP_%i" % (i+1) # Set up array names based on the current model
	
	## Split data into classifier and data 
	
        X = Combined_gene[genelist]
    
        Y = Combined_gene["gene_status"]
	
	## Set up test set and training set
	
        test_size = 0.20
        X_train_train, X_train_valid, y_train_train, y_train_valid = model_selection.train_test_split(X, Y, test_size = test_size)
        eval_set = [(X_train_valid, y_train_valid)]
	
	## Build the model (note colsample_bytree - this limits the number of the genes each model gets, eg. only 30% of all genes each time, leading to a massive reduction in model accuracy if it relies on one gene)
	
        model = XGBClassifier(max_depth = 7, n_estimators = 400, learning_rate=0.1, colsample_bytree = 0.30) ## Columnsampling is here
        model.fit(X_train_train, y_train_train, eval_metric="logloss", early_stopping_rounds=40, eval_set = eval_set, verbose = False) ## Evaluation metric is here
        # Generate weights and import to pandas
        weightscore = model.get_booster().get_score(importance_type = "weight")
        gainscore = model.get_booster().get_score(importance_type = "gain")
        weightframe = importance_frame = pd.DataFrame({importance_array: list(weightscore.values()), 'Feature': list(weightscore.keys())})
        gainframe = importance_frame = pd.DataFrame({importance_array: list(gainscore.values()), 'Feature': list(gainscore.keys())})

        combined_importance = weightframe[importance_array]*gainframe[importance_array]
        combinedframe = pd.concat([gainframe['Feature'],combined_importance], axis=1)

        # Generate SHAP score and import to pandas
        explainer = shap.TreeExplainer(model)
        shap_values = explainer.shap_values(X)
        shap_mean = pd.Series(np.mean(abs(shap_values), axis = 0))
        shap_mean = shap_mean.rename(shap_array)

        pathway_genes = pd.Series(X.columns)
        pathway_genes = pathway_genes.rename("Feature")
        shap_score = pd.concat([shap_mean, pathway_genes], axis = 1)
        # Set up arrays to append total scores
        if i == 0:
            shaprepplot = (pathway + "representative_shap_scores.png")
            shap.summary_plot(shap_values, X, max_display = 25, show = False)
            plt.savefig(shaprepplot, dpi = 300)
            plt.close()
        if i == 0:
            weights_data = weightframe
            gains_data = gainframe
            combined_data = combinedframe
            shap_data = shap_score
        else:
            weights_data = pd.merge(weights_data, weightframe, on='Feature', how="outer")
            gains_data = pd.merge(gains_data, gainframe, on='Feature', how="outer")
            combined_data = pd.merge(combined_data, combinedframe, on='Feature', how="outer")
            shap_data = pd.merge(shap_data, shap_score, on='Feature', how="outer")
            
        ## Calculate overall accuracy of the model for an individual pathway 
        y_pred = model.predict(X_train_valid)
        predictions = [round(value) for value in y_pred]
	# Generate accuray
        accuracy = (accuracy_score(y_train_valid, predictions) * 100) # standard accuracy (percentage correct)
        mcc = (matthews_corrcoef(y_train_valid, predictions) *100) # Matthews Correlation Coefficient (takes into account false positives and false negatives)
        
	# Small function to calculate exact numbers of false pos + neg
        Testset = y_train_valid.as_matrix()
        zipped = np.dstack((Testset, predictions))    
        for item in zipped:
             falsepos = 0
             falseneg = 0
             match = 0     
             for place in item:
                 totalcount = len(item)
                 if (place[0] == 0) & (place[1] == 1):
                     falsepos += 1
                 elif (place[0] == 1) & (place[1] == 0):
                     falseneg += 1
                 elif place[0] == place[1]:
                     match += 1
        for item in zipped:
            negcount = 0
            poscount = 0
            for place in item:
                if place[0] == 0:
                    negcount += 1
                if place[0] == 1:
                    poscount += 1

                        
        if not negcount == 0:
            falseposp = (falsepos/negcount) * 100
        else:
            falseposp = "Invalid"
        if not poscount == 0:
            falsenegp = (falseneg/poscount) * 100
        else:
            falsenegp = "Invalid"
            
        totalp = (match/totalcount) * 100
	# Append accuracies to dataframe
        accuracy_frame = pd.Series([accuracy, falseposp, falsenegp, mcc], name = i, index = ["Accuracy", "Falsepos", "Falseneg", "MCC_score"])
        if i == 0:
            accuracy_base = accuracy_frame.to_frame()
        else:
            accuracy_base = accuracy_base.join(accuracy_frame, how = "outer")
        ## Loop ends - now have dataframes containing scores for each model, for each pathway

    ## Do some general curation of the dataframes
    
    # Clean dataframes (add zeros in place of NA)
    weights_data = weights_data.fillna(0)
    gains_data = gains_data.fillna(0)
    combined_data = combined_data.fillna(0)
    shap_data = shap_data.fillna(0)
    # Rename indexes
    weights_data = weights_data.set_index('Feature')
    gains_data = gains_data.set_index('Feature')
    combined_data = combined_data.set_index('Feature')
    shap_data = shap_data.set_index('Feature')

    # Calculate averages and standard deviation
    weights_data["Average"] = weights_data.mean(axis =1)
    weights_data["STD"] = weights_data.std(axis = 1)

    gains_data["Average"] = gains_data.mean(axis =1)
    gains_data["STD"] = gains_data.std(axis = 1)

    combined_data["Average"] = combined_data.mean(axis =1)
    combined_data["STD"] = combined_data.std(axis = 1)

    shap_data["Average"] = shap_data.mean(axis =1)
    shap_data["STD"] = shap_data.std(axis =1)

    # Do some sorting
    weights_data.sort_values(by = "Average", inplace = True, ascending = False)
    gains_data.sort_values(by = "Average", inplace = True, ascending = False)
    combined_data.sort_values(by = "Average", inplace = True, ascending = False)
    shap_data.sort_values(by = "Average", inplace = True, ascending = False)

    ## Plot feature weights/gains
    
    # Shorten data to avoid huge plots
    if len(combined_data.index) > 25: 
        combined_data = combined_data.head(25)
        weights_data = weights_data.head(25)
        gains_data = gains_data.head(25)
        shap_data = shap_data.head(25) 
    	
    # Generate plot names
    weightname = (pathway + "featureweight_average.png")
    gainname = (pathway + "featuregain_average.png")
    combinedname = (pathway + "featurecombined_average.png")
    shapname = (pathway + "shapcombined_average.png")
    
    # Plot weights
    palette = sns.light_palette("#43A9DB", n_colors=len(combined_data.index), reverse = True)
    sns.barplot(y = weights_data.index, x=weights_data.Average, palette= palette, edgecolor = "k", xerr = weights_data.STD)
    sns.despine()
    plt.tight_layout()
    plt.savefig(weightname, dpi = 300)
    plt.close()

    # Plot gains
    sns.barplot(y = gains_data.index, x=gains_data.Average, palette= palette, edgecolor = "k", xerr = gains_data.STD)
    sns.despine()
    plt.tight_layout()
    plt.savefig(gainname, dpi = 300)
    plt.close()

    # Plot weights * gains
    sns.barplot(y = combined_data.index, x=combined_data.Average, palette= palette, edgecolor = "k", xerr = combined_data.STD)
    sns.despine()
    plt.tight_layout()
    plt.savefig(combinedname, dpi = 300)
    plt.close()

    # Plot SHAP accuracy
    sns.barplot(y = shap_data.index, x=shap_data.Average, palette= palette, edgecolor = "k", xerr = shap_data.STD)
    sns.despine()
    plt.tight_layout()
    plt.savefig(shapname, dpi = 300)
    plt.close()
    
    # Calculate averages for all models
    accuracy_base["Average"] = accuracy_base.mean(axis =1)
    accuracy_base["STD"] = accuracy_base.std(axis = 1)


    # Generate and manipulate stats array
    pathwayaverage = pathway + "_average"
    pathwaystd = pathway + "_std"
    averagevalues = accuracy_base["Average"].to_frame().transpose()
    averagevalues.columns = ["Av_acc", "Av_Falsepos", "Av_Falseneg", "Av_MCC"]
    stdvalues = accuracy_base["STD"].to_frame().transpose()
    stdvalues.columns = ["Std_acc", "Std_Falsepos", "Std_Falseneg", "Std_MCC"]
    averagevalues.reset_index(drop = True, inplace = True)
    stdvalues.reset_index(drop = True, inplace = True)
    combined_analysis = pd.concat([averagevalues, stdvalues], axis = 1)
    combined_analysis["Pathway"] = pathway
    pathwaystats.append(combined_analysis)


    ## Save all data to files
    
    weightframename = (pathway + "featureweight_values.csv")
    gainframename = (pathway + "featuregain_values.csv")
    combinedframename = (pathway + "featurecombined_values.csv")
    accuracy_framename = (pathway + "accuracy_scores.csv")
    shap_framename = (pathway + "shap_scores.csv")

    pd.DataFrame.to_csv(weights_data, weightframename)
    pd.DataFrame.to_csv(gains_data, gainframename)
    pd.DataFrame.to_csv(combined_data, combinedframename)
    pd.DataFrame.to_csv(accuracy_base, accuracy_framename) 
    pd.DataFrame.to_csv(shap_data, shap_framename) 
    os.chdir("../")
    
## Save all pathways to a single array

outputarray = pd.concat(pathwaystats, axis = 0)
outputarray.sort_values(by = "Av_acc", inplace = True, ascending = False)

pd.DataFrame.to_csv(outputarray, "All_pathway_stats.csv")

## Plot accuracy and MCC scores for the whole 
sns.barplot(y = outputarray.Pathway, x=outputarray.Av_acc, palette= palette, edgecolor = "k", xerr = outputarray.Std_acc)
plt.tight_layout()
plt.xlim([0,100])
sns.despine()
plt.savefig("Pathway_accuracy.png", dpi = 300)
plt.close()

sns.barplot(y = outputarray.Pathway, x=outputarray.Av_MCC, palette= palette, edgecolor = "k", xerr = outputarray.Std_MCC)
plt.tight_layout()
plt.xlim([0,100])
sns.despine()
plt.savefig("Pathway_MCC.png", dpi = 300)
plt.close()
