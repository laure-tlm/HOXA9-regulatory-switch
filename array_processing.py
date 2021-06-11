## Script for processing TCGA data - data is downloaded from firebrowse (http://firebrowse.org/). Two files are required - one ending in RSEM_genes (contains raw count and pre-TPM data), and RSEM_genes_normalized (contains RSEM upper quartile normalised data). This script processes those files for easier analysis.
import pandas as pd
import glob
import argparse, argcomplete
import os, shutil

### Do all the parsing stuff ###
parser = argparse.ArgumentParser()
parser.add_argument("-mode", dest= "mode", choices=["test", "live"], default= "live")
argcomplete.autocomplete(parser)
args = parser.parse_args()

if args.mode == "test":
	print("#### Script is running in test mode - it will overwrite all directories ####")



# Generate Data Directories
if args.mode == "test":
	for foldername in ("../TPM_data", "../RSEM_upperquartile_norm_data", "../Raw_counts_data", "../Projects"):
		if os.path.exists(foldername):
			shutil.rmtree(foldername)
		os.mkdir(foldername)
	for name in glob.glob("./**/*genes_normalized__data.data.txt", recursive = True):
		datasetname = (((os.path.basename(name)).split(".")[0]))
		os.mkdir("../Projects/" + datasetname)

else:
	for foldername in ("../TPM_data", "../RSEM_upperquartile_norm_data", "../Raw_counts_data", "../Projects"):
		if not os.path.exists(foldername):
			os.mkdir(foldername)
	for name in glob.glob("./**/*genes_normalized__data.data.txt", recursive = True):
		datasetname = (((os.path.basename(name)).split(".")[0]))
		if not os.path.exists("../Projects/" + datasetname):
			os.mkdir("../Projects/" + datasetname)


# Search directory for data files
for name in glob.glob("./**/*.data.txt", recursive = True):
	# Take the first dataset
	if name.endswith("genes__data.data.txt"):
		
		# Generate name of TCGA project
		datasetname = (((os.path.basename(name)).split(".")[0]))
		print("Processing " + datasetname + " dataset")
		data1 = pd.read_csv(name, sep="\t", low_memory = True)

		# Remove useless columns and split into raw counts and normalised counts
		data1_dropped = data1.drop((data1.columns[data1.iloc[0] == "transcript_id"]), axis =1)
		data1_rawcounts = data1_dropped.drop((data1_dropped.columns[data1_dropped.iloc[0] == "scaled_estimate"]), axis =1)
		data1_scaledestimate = data1_dropped.drop((data1_dropped.columns[data1_dropped.iloc[0] == "raw_count"]), axis =1)

		# Clean up the tables a bit
		data1_rawcounts_final = data1_rawcounts.set_index("Hybridization REF")
		data1_rawcounts_final = data1_rawcounts_final.drop("gene_id")
		data1_scaledestimate_final = data1_scaledestimate.set_index("Hybridization REF")
		data1_scaledestimate_final = data1_scaledestimate_final.drop("gene_id")

		# Convert scaled estimate to TPM
		print("Performing TPM transformation")
		data1_float = data1_scaledestimate_final.apply(pd.to_numeric)
		data1_TPM = data1_float.mul(1000000)

		# Save files
		print("Saving " + datasetname + " data to csv")
		rawcounts_filename = (datasetname + "_raw_counts.csv")
		TPM_filename = (datasetname + "_TPM.csv")
		data1_rawcounts_final.to_csv("../Raw_counts_data/" + rawcounts_filename)
		data1_rawcounts_final.to_csv("../Projects/" + datasetname + "/" + rawcounts_filename)
		data1_TPM.to_csv("../TPM_data/" + TPM_filename)
		data1_TPM.to_csv("../Projects/" + datasetname + "/" + TPM_filename)

	elif name.endswith("genes_normalized__data.data.txt"):

		# Generate name of TCGA project
		datasetname = (((os.path.basename(name)).split(".")[0]))
		print("Processing " + datasetname + " normalised dataset")
		data2 = pd.read_csv(name, sep="\t", low_memory = True)

		# Clean up tables
		data2_final = data2.set_index("Hybridization REF")
		data2_final = data2_final.drop("gene_id")

		# Save files
		print("Saving " + datasetname + " data to csv")
		RSEM_upperquartile_filename = (datasetname + "_RSEM_upperquartile_norm.csv")
		data2_final.to_csv("../RSEM_upperquartile_norm_data/" + RSEM_upperquartile_filename)
		data2_final.to_csv("../Projects/" + datasetname + "/" + RSEM_upperquartile_filename)

