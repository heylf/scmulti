
# In[]
print("[START]")

import argparse
import muon as mu
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import json

#########################################
###### GLOBAL VARS and DIRECTORIES ######
#########################################

tool_description = """
Runnin QC
"""

# parse command line arguments
parser = argparse.ArgumentParser(description=tool_description, formatter_class=argparse.RawDescriptionHelpFormatter)

# version
parser.add_argument("-o", "--out", dest="out", metavar='str', required=True, help="Single cell files")
parser.add_argument("-v", "--version", action="version", version="%(prog)s 0.1.0")
parser.add_argument("-s", "--samples", dest="samples", metavar='str', required=True, help="Sample names")
parser.add_argument("-f", "--files", dest="files", metavar='str', required=True, help="Output dir")

args = vars(parser.parse_args())

#######################
###### LOAD DATA ######
#######################
# This will get us the data from cellranger-arc
# It also generates a table for the necessary information (arc_aggr_df.to_csv)
print("[TASK] Load data")

atacs = []
rnas = []
samples = args['samples'].split(",")

num_samples = 0

for file in args['files'].split(","):
    print(file)

    cdata = mu.read_10x_h5(file)

    cdata["rna"].var_names_make_unique()
    cdata["atac"].var_names_make_unique()
    
    atacs.append(cdata["atac"])
    rnas.append(cdata["rna"])

sample_name_sorted_indices = np.argsort(np.array(samples))

# Get RNA counts
rna = ad.concat([rnas[i] for i in sample_name_sorted_indices], 
                label="sample", keys=[samples[i] for i in sample_name_sorted_indices], index_unique="_", merge="same")

# Get ATAC counts
atac = ad.concat([atacs[i] for i in sample_name_sorted_indices], 
                 label="sample", keys=[samples[i] for i in sample_name_sorted_indices], index_unique="_", merge="same")

########################
###### PREPARE QC ######
########################
print("[TASK] Prepare QC")

# Perform QC
rna.var["mt"] = rna.var_names.str.startswith("MT-") # this will add mitochondrial QC to the general QC
rna.var["ribo"] = rna.var_names.str.startswith("RPS") | rna.var_names.str.startswith("RPL") # this will add ribosomal RNA QC
sc.pp.calculate_qc_metrics(rna, qc_vars=["mt", "ribo"], inplace=True)

# Add QC levels for samples
sample_series = rna.obs.groupby("sample")

###########################
###### MULTIQC PLOTS ######
###########################
print("[TASK] Generate Multiqc plots")

d = rna.obs.groupby("sample").size().sort_values(ascending=False).rename("ncells").reset_index().assign(sample=lambda df: pd.Categorical(df["sample"].to_list(), df["sample"].to_list(), ordered=True))

data = {}
for index, row in d.iterrows():
    data[row['sample']] = {'cells': row['ncells']}

# Your JSON data
json_data = {
    "id": "custom_data_lineplot",
    "section_name": "Custom JSON File",
    "description": "This plot is a self-contained JSON file.",
    "plot_type": "bargraph",
    "pconfig": {
        "id": "custom_data_linegraph",
        "title": "Output from my JSON file",
        "ylab": "Number of things",
        "xDecimals": False
    },
    "data": data 
}

# Writing JSON data to a file
with open(args['out'], 'w') as json_file:
    json.dump(json_data, json_file, indent=4)

print("[FINISH]")
# %%
