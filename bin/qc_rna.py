#!/usr/bin/env python

# In[]

import argparse
import muon as mu
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import json
import plotly.express as px
import plotly.graph_objects as go
import random

# pip install nbformat

# In[]
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

# In[]
print("[START]")

# Seeds (!!! DO NOT CHANGE THIS SEED !!!)
seed = 123
random.seed(seed)
print(f"[NOTE] seed {seed}")


# Calculate confidence intervale function
def get_ci_df(count_series, mode, rna):
    min_num_cells = np.min(count_series.size())

    # Get number of groups
    num_group = len(rna.obs.groupby(mode).size())

    num_replicates = 1000

    mean_array = np.empty([0])
    labels = np.empty([0])

    for key, item in count_series:
        labels = np.append(labels, np.array([key]*num_replicates))

        tmp_mean_array = np.empty([num_replicates])
        counts = np.array(count_series.get_group(key).tolist())

        for i in range(0, num_replicates):
            tmp_mean_array[i] = np.mean(random.choices(counts, k=min_num_cells))

        mean_array = np.append(mean_array, tmp_mean_array)
    
    d = {'group': labels, 'values': mean_array}
    d = pd.DataFrame(data=d, index=[x for x in range(0, num_group*num_replicates)])
    return(d)

#######################
###### LOAD DATA ######
#######################
# This will get us the data from cellranger-arc
# It also generates a table for the necessary information (arc_aggr_df.to_csv)
print("[TASK] Load data")

atacs = []
rnas = []
# samples = args['samples'].split(",")
samples = 'test_ARC'

num_samples = 0

# for file in args['files'].split(","):
for file in ['/home/florian/Documents/tmp_data_folder/output/cellrangerarc/count/test_scARC/outs/filtered_feature_bc_matrix.h5']:
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

#In[]
# TODO add thresholds for UMI etc.

###########################
###### MULTIQC PLOTS ######
###########################
print("[TASK] Generate Multiqc plots")

d = rna.obs.groupby("sample").size().sort_values(ascending=False).rename("ncells").reset_index().assign(sample=lambda df: pd.Categorical(df["sample"].to_list(), df["sample"].to_list(), ordered=True))


figures = []

# ----------------------------------------------------------------------------------------------------------------------
# Cell Counts ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
fig = px.bar(d, x="sample", y="ncells", text="ncells", barmode="stack")
fig.update_layout(
    yaxis_title="Number of cells",
    xaxis_title="Samples",
    title="Total number of cells per sample",
    hovermode=False  # Deactivate hover info
)
fig.update_layout(xaxis_tickangle=-90)
figures.append(fig)


# ----------------------------------------------------------------------------------------------------------------------
# UMI/Gene/MT/RB -------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

for level in ["total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo"]:

    fig = px.histogram(rna.obs, x=level, nbins=100,
                    labels={level: f"{level}"},
                    width=800, height=800)
    fig.update_layout(title=f"Distribution of {level} for all samples",
                    showlegend=True)
    figures.append(fig)


    # TODO make more beautiful
    fig = px.histogram(rna.obs, x=level, nbins=100,
                    facet_col="sample", facet_col_wrap=4,
                    log_x=True, labels={level: f"log10 {level}"},
                    width=800, height=800)
    fig.update_layout(title=f"Distribution of log10 {level}",
                    showlegend=True)
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    figures.append(fig)


    fig = px.histogram(rna.obs, x=level, nbins=100,
                    facet_col="sample", facet_col_wrap=4,
                    labels={level: f"{level}"},
                    width=800, height=800)
    fig.update_layout(title=f"Distribution of {level}",
                    showlegend=True)
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    figures.append(fig)

    # Get the data for the boxplot
    ci_df = get_ci_df(rna.obs.groupby("sample")[level], "sample", rna)

    # Create the boxplot trace
    boxplot_trace = go.Box(
        x=ci_df["group"],
        y=ci_df["values"],
        name="Means of total counts"
    )

    # Create the layout
    layout = go.Layout(
        xaxis=dict(
            tickangle=-45,
            tickfont=dict(size=10),
            title="sample"
        ),
        yaxis=dict(
            title=f"Means of {level}"
        ),
        height=500,
        width=800
    )

    # Create the figure
    fig = go.Figure(data=[boxplot_trace], layout=layout)
    figures.append(fig)

# ----------------------------------------------------------------------------------------------------------------------
# Others -------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

THRESH_GENE_FILTER = 100

num_genes = len([1 for i in rna.var.n_cells_by_counts > min(THRESH_GENE_FILTER, np.ceil(0.01 * rna.n_obs)) if i == False])

fig = px.histogram(rna.var, x="n_cells_by_counts", nbins=100,
                   labels={"n_cells_by_counts": "log10 number of cells expressing > 0 (n reads)"},
                   log_x=True,
                   title=f"Distribution of the number of cells per gene. From total number of {rna.n_vars} genes, {num_genes} genes have < {THRESH_GENE_FILTER} cells.")
fig.add_vline(x=THRESH_GENE_FILTER)
figures.append(fig)


THRESH_UMI_N_GENES = 500
THRESH_MT = 20

fig = go.Figure()

fig.add_trace(go.Scatter(
    x=rna.obs["total_counts"],
    y=rna.obs["n_genes_by_counts"],
    mode="markers",
    marker=dict(
        size=rna.obs["pct_counts_ribo"],
        sizemode="diameter",
        sizeref=0.1,
        sizemin=1,
        color=rna.obs["pct_counts_mt"],
        colorscale="Viridis",
        cmin=0,
        cmax=100,
        colorbar=dict(
            title="pct_counts_mt",
            titleside="right"
        )
    )
))

fig.update_layout(
    title="Scatter plot of total_counts vs n_genes_by_counts vs pct_counts_mt vs pct_counts_ribo",
    xaxis=dict(
        type="log",
        title="log10 total_counts"
    ),
    yaxis=dict(
        type="log",
        title="log10 n_genes_by_counts"
    ),
    showlegend=False
)

if fig.layout.yaxis.range is not None:
    fig.add_shape(
        type="line",
        x0=THRESH_UMI_N_GENES,
        y0=fig.layout.yaxis.range[0],
        x1=THRESH_UMI_N_GENES,
        y1=fig.layout.yaxis.range[1],
        line=dict(
            color="black",
            width=1,
            dash="dash"
        )
    )

if fig.layout.yaxis.range is not None:
    fig.add_shape(
        type="line",
        x0=fig.layout.xaxis.range[0],
        y0=THRESH_UMI_N_GENES,
        x1=fig.layout.xaxis.range[1],
        y1=THRESH_UMI_N_GENES,
        line=dict(
            color="black",
            width=1,
            dash="dash"
        )
    )
figures.append(fig)


# Plot the highest expressed genes before filtering
import plotly.graph_objects as go

nGENES = 30

idx = np.argsort(rna.var["n_cells_by_counts"])[-nGENES:]

expression = rna.X.todense()[:, idx]
genes = rna.var_names[idx]

df = pd.DataFrame(expression, columns=genes)

fig = go.Figure()

for gene in genes:
    fig.add_trace(go.Box(y=df[gene], name=gene))

fig.update_layout(
    title=f"Boxplot of {nGENES} highest expressed genes",
    xaxis=dict(title="Genes"),
    yaxis=dict(title="Expression"),
)

figures.append(fig)

# ----------------------------------------------------------------------------------------------------------------------
# Generate HTML --------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

with open('rna_qc_mqc.html', 'w') as f:
    for fig in figures:
        f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))
    
print("[FINISH]")
# %%
