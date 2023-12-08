#!/usr/bin/env python


import argparse
import muon as mu
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import random

# pip install nbformat
# pip install pysam
# pip install plotly

########################################################################################################################
###### GLOBAL VARS and DIRECTORIES #####################################################################################
########################################################################################################################

tool_description = """
Runnin QC
"""

# parse command line arguments
parser = argparse.ArgumentParser(description=tool_description, formatter_class=argparse.RawDescriptionHelpFormatter)

# version
parser.add_argument("-v", "--version", action="version", version="%(prog)s 0.1.0")

# Required arguments
parser.add_argument("-s", "--samples", dest="samples", metavar='str', required=True, help="Sample names")
parser.add_argument("-f", "--files", dest="files", metavar='str', required=True, help="h5 files")

# Optional arguments
parser.add_argument("-o", "--out", dest="out", metavar='str', required=False, help="Single cell files")
parser.add_argument("--demux", dest="demux", metavar='str', required=False, help="Demultiplexed files")
parser.add_argument("--thresh_gene", dest="thresh_gene", metavar='str', required=False, 
                    help="Threshold for gene cell count filter", default=100)
parser.add_argument("--thresh_umi_genecount", dest="thresh_umi_genecount", metavar='str', required=False, 
                    help="Threshold for UMI and gene counr filter", default=500)
parser.add_argument("--thresh_mt", dest="thresh_mt", metavar='str', required=False, 
                    help="Threshold for mitochondrial filter", default=20)
parser.add_argument("--thresh_rb", dest="thresh_rb", metavar='str', required=False, 
                    help="Threshold for ribsomal coverage filter", default=10)

args = vars(parser.parse_args())

THRESH_GENE_FILTER = args["thresh_gene"]
THRESH_UMI_N_GENES = args["thresh_umi_genecount"]
THRESH_MT = args["thresh_mt"]
THRESH_RB = args["thresh_rb"]

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

########################################################################################################################
###### LOAD DATA #######################################################################################################
########################################################################################################################
# This will get us the data from cellranger-arc
# It also generates a table for the necessary information (arc_aggr_df.to_csv)
print("[TASK] Load data")

rnas = []

samples = args['samples'].split(",")

for file in args['files'].split(","):
    print(file)

    cdata = mu.read_10x_h5(file)

    cdata["rna"].var_names_make_unique()
    
    rnas.append(cdata["rna"])

sample_name_sorted_indices = np.argsort(np.array(samples))

# Get RNA counts
rna = ad.concat([rnas[i] for i in sample_name_sorted_indices], 
                label="sample", keys=[samples[i] for i in sample_name_sorted_indices], index_unique="_", merge="same")

########################################################################################################################
###### PREPARE QC ######################################################################################################
########################################################################################################################
print("[TASK] Prepare QC")

# Perform QC
rna.var["mt"] = rna.var_names.str.startswith("MT-") # this will add mitochondrial QC to the general QC
rna.var["ribo"] = rna.var_names.str.startswith("RPS") | rna.var_names.str.startswith("RPL") # this will add ribosomal RNA QC
sc.pp.calculate_qc_metrics(rna, qc_vars=["mt", "ribo"], inplace=True)

########################################################################################################################
###### QC PLOTS RNA ####################################################################################################
########################################################################################################################
print("[TASK] Generate QC plots RNA")

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
    fig.update_layout(title=f"Total distribution of {level} for all samples",
                    showlegend=True)
    figures.append(fig)

    if ( level not in ["pct_counts_mt", "pct_counts_ribo"] ):

        rna.obs[f'log_{level}'] = np.log10(rna.obs[level] + 1)

        # Create a histogram trace for each sample
        traces = []
        for sample in samples:
            sample_data = rna.obs[rna.obs['sample'] == sample]
            trace = go.Histogram(x=sample_data[f'log_{level}'], nbinsx=100, name=sample)
            traces.append(trace)

        # Create the merged plot
        fig = go.Figure(data=traces)

        # Update layout
        fig.update_layout(
            title=f"Distribution of log10 {level} for all samples",
            showlegend=True,
            height=800,
            margin=dict(t=150)  # Add more space between title and plot
        )

        fig.update_layout(
            dict(
                updatemenus=[
                    dict(
                        type="buttons",
                        direction="left",
                        buttons=list([
                            dict(
                                args=["visible", "legendonly"],
                                label="Deselect All",
                                method="restyle"
                            ),
                            dict(
                                args=["visible", True],
                                label="Select All",
                                method="restyle"
                            )
                        ]),
                        pad={"r": 10, "t": 10},
                        showactive=False,
                        x=1,
                        xanchor="right",
                        y=1.1,
                        yanchor="top"
                    ),
                ]
            )
        )

        fig.update_layout(
            xaxis_title=f"log_{level}",
            yaxis_title="Frequency"
        )

        figures.append(fig)

    # Create a histogram trace for each sample
    traces = []
    for sample in samples:
        sample_data = rna.obs[rna.obs['sample'] == sample]
        trace = go.Histogram(x=sample_data[level], nbinsx=100, name=sample)
        traces.append(trace)

    # Create the merged plot
    fig = go.Figure(data=traces)

    # Update layout
    fig.update_layout(
        title=f"Distribution of {level} for all samples",
        showlegend=True,
        height=800,
        margin=dict(t=150)  # Add more space between title and plot
    )

    fig.update_layout(
        dict(
            updatemenus=[
                dict(
                    type="buttons",
                    direction="left",
                    buttons=list([
                        dict(
                            args=["visible", "legendonly"],
                            label="Deselect All",
                            method="restyle"
                        ),
                        dict(
                            args=["visible", True],
                            label="Select All",
                            method="restyle"
                        )
                    ]),
                    pad={"r": 10, "t": 10},
                    showactive=False,
                    x=1,
                    xanchor="right",
                    y=1.1,
                    yanchor="top"
                ),
            ]
        )
    )

    fig.update_layout(
        xaxis_title=f"{level}",
        yaxis_title="Frequency"
    )

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
            title="Sample"
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

num_genes = len([1 for i in rna.var.n_cells_by_counts > min(THRESH_GENE_FILTER, np.ceil(0.01 * rna.n_obs)) if i == False])

rna.var['log_n_cells_by_counts'] = np.log10(rna.var['n_cells_by_counts'] + 1)

fig = px.histogram(rna.var, x="log_n_cells_by_counts", nbins=100,
                   labels={"log_n_cells_by_counts": "log10 number of cells expressing > 0 (n reads)"},
                   title=f"Distribution of the number of cells per gene. From total number of {rna.n_vars} genes, {num_genes} genes have < {THRESH_GENE_FILTER} cells.")
fig.add_vline(x=np.log10(THRESH_GENE_FILTER))
figures.append(fig)

fig = go.Figure()

fig.add_trace(go.Scatter(
    x=rna.obs["total_counts"],
    y=rna.obs["n_genes_by_counts"],
    mode="markers",
    marker=dict(
        size=rna.obs["pct_counts_ribo"],
        sizemode="diameter",
        sizeref=5,
        sizemin=0.1,
        color=rna.obs["pct_counts_mt"],
        colorscale="RdYlBu_r",
        cmin=0,
        cmax=100,
        colorbar=dict(
            title="pct_counts_mt",
            titleside="right"
        ),
        showscale=True,  # Add this line to show the size scale legend
        colorbar_title="pct_counts_mt"
    ),
    text=[
        f"total_counts: {count}<br>n_genes_by_counts: {genes}<br>pct_counts_ribo: {ribo}<br>pct_counts_mt: {mt}"
        for count, genes, ribo, mt in zip(
            rna.obs["total_counts"],
            rna.obs["n_genes_by_counts"],
            np.round(rna.obs["pct_counts_ribo"]),
            np.round(rna.obs["pct_counts_mt"])
        )
    ]
))

fig.add_shape(
    type="line",
    x0=THRESH_UMI_N_GENES,
    x1=THRESH_UMI_N_GENES,
    y1=max(rna.obs["n_genes_by_counts"]),
    line=dict(
        color="black",
        width=1
    )
)

fig.add_shape(
    type="line",
    y0=THRESH_UMI_N_GENES,
    y1=THRESH_UMI_N_GENES,
    x1=max(rna.obs["total_counts"]),
    line=dict(
        color="black",
        width=1
    )
)

fig.update_layout(
    title="Scatter plot of total_counts vs n_genes_by_counts vs pct_counts_mt vs pct_counts_ribo. <br>\
        Size of the diameter corresponds to the pct_counts_ribo.",
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

figures.append(fig)

# Plot the highest expressed genes before filtering
import plotly.graph_objects as go

nGENES = 30

idx = np.argsort(rna.var["n_cells_by_counts"])[-nGENES:]

expression = rna.X.todense()[:, idx]
genes = rna.var_names[idx]

df = pd.DataFrame(expression, columns=genes)
df_perct = df.div(np.array(rna.obs['total_counts']), axis=0)

# Calculate the median for each gene
gene_means = df_perct.mean()

# Sort the genes based on their median values
sorted_genes = gene_means.sort_values().index[::-1]

fig = go.Figure()

for gene in sorted_genes:
    fig.add_trace(go.Box(y=df_perct[gene], name=gene))

fig.update_layout(
    title=f"Boxplot of {nGENES} highest expressed genes",
    xaxis=dict(title="Genes"),
    yaxis=dict(title="% of total counts")
)

figures.append(fig)

# ----------------------------------------------------------------------------------------------------------------------
# Generate HTML --------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

with open(f'rna_qc_sample_mqc.html', 'w') as f:
    for fig in figures:
        f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))


########################################################################################################################
###### QC PLOTS RNA DEMULIPLEXED #######################################################################################
########################################################################################################################

if ( args['demux'] ):

    demux_files = args['demux'].split(",")

    # Do Demultiplexing
    assignments = []
    for i, file in enumerate(demux_files):
        ass = pd.read_table(file)[["cell", "donor_id"]]
        ass.cell += f"_{samples[i]}"
        ass.set_index("cell", inplace=True)
        assignments.append(ass)
    assignments = pd.concat(assignments, axis=0)

    rna.obs = rna.obs.join(assignments, how="left")

    # Get number of donors
    ndonor = len(set(rna.obs["donor_id"].tolist()))

    # Generate Colors for organoids
    r = lambda: random.randint(0,255)
    DONOR_COLORS = ['#%02X%02X%02X' % (r(),r(),r()) for i in range(0,ndonor)]
    DONOR_COLORS[sorted(set(rna.obs["donor_id"].tolist())).index("unassigned")] = "Black"
    DONOR_COLORS[sorted(set(rna.obs["donor_id"].tolist())).index("doublet")] = "Red"


    figures = []

    donor_series = rna.obs.groupby("sample")["donor_id"]
    num_donors = []
    num_samples = len(set(rna.obs["sample"].tolist()))

    fractions = [0.0, 0.0, 0.0] *  num_samples
    fractions_cat = ["unassigned", "doublet", "assigned"] * num_samples

    i = 0
    n = 0
    for i, item in enumerate(donor_series):

        sample = item[0]
            
        sample_donor_df = donor_series.get_group(sample).value_counts()

        total_cell = sum(sample_donor_df)

        fractions[i + n] = sample_donor_df["unassigned"]/total_cell
        fractions[i + (n+1)] = sample_donor_df["doublet"]/total_cell

        sample_donor_df = sample_donor_df.drop(["unassigned", "doublet"], axis=0)
        num_donors.append(len(set(sample_donor_df.index.values)))

        fractions[i + (n+2)] = sum(sample_donor_df)/total_cell

        n = i + 2
        i += 1

    d = {'sample': samples, 'num_donors': num_donors}
    d = pd.DataFrame(data=d, index=samples)

    # Plot
    fig = px.bar(d, x='sample', y='num_donors')
    fig.update_layout(
        xaxis=dict(title="Sample"),
        yaxis=dict(title="Number of Donors"),
        title="Number of donors per sample identified by the demultiplexing tool"
    )
    figures.append(fig)

    d = {'sample': [s for s in samples for _ in (0, 1, 2)], 'category_cells': fractions_cat, 
        'frac_cells': fractions}
    d = pd.DataFrame(data=d, index=[x for x in range(0, len(samples)*3)])

    fig = px.bar(d, x='sample', y='frac_cells', color='category_cells')
    fig.update_layout(
        xaxis=dict(title="Sample"),
        yaxis=dict(title="Fraction of cells"),
        title="Fraction of cells per sample assigned by the demultiplexing tool"
    )
    figures.append(fig)

    # Plot
    d = rna.obs.groupby("donor_id").size().sort_values(ascending=False).rename("ncells").reset_index()\
        .assign(donor=lambda df: pd.Categorical(df["donor_id"].to_list(), df["donor_id"].to_list(), ordered=True))
        
    # ------------------------------------------------------------------------------------------------------------------
    # Cell Counts ------------------------------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------------------------------
    fig = px.bar(d, x="donor_id", y="ncells", text="ncells", barmode="stack")
    fig.update_layout(
        yaxis_title="Number of cells",
        xaxis_title="Samples",
        title="Total number of cells per sample",
        hovermode=False  # Deactivate hover info
    )
    fig.update_layout(xaxis_tickangle=-90)
    figures.append(fig)     

    # ------------------------------------------------------------------------------------------------------------------
    # UMI/Gene/MT/RB ---------------------------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------------------------------

    # Get unique donors
    unique_donors = rna.obs['donor_id'].unique()

    for level in ["total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo"]:
    
        # Create a histogram trace for each donor
        traces = []
        for donor in unique_donors:
            donor_data = rna.obs[rna.obs['donor_id'] == donor]
            trace = go.Histogram(x=donor_data[level], nbinsx=100, name=donor)
            traces.append(trace)

        # Create the merged plot
        fig = go.Figure(data=traces)

        # Update layout
        fig.update_layout(
            title=f"Distribution of {level} for all donors",
            showlegend=True,
            height=800,
            margin=dict(t=150)  # Add more space between title and plot
        )

        fig.update_layout(
            dict(
                updatemenus=[
                    dict(
                        type="buttons",
                        direction="left",
                        buttons=list([
                            dict(
                                args=["visible", "legendonly"],
                                label="Deselect All",
                                method="restyle"
                            ),
                            dict(
                                args=["visible", True],
                                label="Select All",
                                method="restyle"
                            )
                        ]),
                        pad={"r": 10, "t": 10},
                        showactive=False,
                        x=1,
                        xanchor="right",
                        y=1.1,
                        yanchor="top"
                    ),
                ]
            )
        )

        fig.update_layout(
            xaxis_title=f"{level}",
            yaxis_title="Frequency"
        )

        figures.append(fig)

        if ( level not in ["pct_counts_mt", "pct_counts_ribo"] ):

            # Create a histogram trace for each donor
            traces = []
            for donor in unique_donors:
                donor_data = rna.obs[rna.obs['donor_id'] == donor]
                trace = go.Histogram(x=donor_data[f'log_{level}'], nbinsx=100, 
                                    name=donor)
                traces.append(trace)

            # Create the merged plot
            fig = go.Figure(data=traces)

            # Update layout
            fig.update_layout(
                title=f"Distribution of log10 {level} for all donors",
                showlegend=True,
                height=800,
                margin=dict(t=150)  # Add more space between title and plot
            )

            fig.update_layout(
                dict(
                    updatemenus=[
                        dict(
                            type="buttons",
                            direction="left",
                            buttons=list([
                                dict(
                                    args=["visible", "legendonly"],
                                    label="Deselect All",
                                    method="restyle"
                                ),
                                dict(
                                    args=["visible", True],
                                    label="Select All",
                                    method="restyle"
                                )
                            ]),
                            pad={"r": 10, "t": 10},
                            showactive=False,
                            x=1,
                            xanchor="right",
                            y=1.1,
                            yanchor="top"
                        ),
                    ]
                )
            )

            fig.update_layout(
                xaxis_title=f"Log10 {level}",
                yaxis_title="Frequency"
            )

            figures.append(fig)

        # Get the data for the boxplot
        ci_df = get_ci_df(rna.obs.groupby("donor_id")[level], "donor_id", rna)

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
                title="Donor"
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

    # ------------------------------------------------------------------------------------------------------------------
    # Sample Donor Level -----------------------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------------------------------

    for level in ["log_total_counts", "log_n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo"]:

        from plotly.subplots import make_subplots

        # Create a subplot with multiple histograms
        fig = make_subplots(rows=1, cols=len(samples), subplot_titles=samples)

        # Create a histogram trace for each donor within each sample subplot
        for i, sample in enumerate(samples):
            sample_data = rna.obs[rna.obs['sample'] == sample]
            
            traces = []
            for donor in unique_donors:
                donor_data = sample_data[sample_data['donor_id'] == donor]
                trace = go.Histogram(x=donor_data[level], nbinsx=100, name=f'{sample} - {donor}')
                traces.append(trace)
            
            # Add traces to the subplot
            for trace in traces:
                fig.add_trace(trace, row=1, col=i + 1)

        # Update layout
        fig.update_layout(
            title=f"Distribution of {level} for all samples and the corresponding donors",
            showlegend=True,
            height=800,
            margin=dict(t=150),  # Add more space between title and plot
        )

        fig.update_layout(
            dict(
                updatemenus=[
                    dict(
                        type="buttons",
                        direction="left",
                        buttons=list([
                            dict(
                                args=["visible", "legendonly"],
                                label="Deselect All",
                                method="restyle"
                            ),
                            dict(
                                args=["visible", True],
                                label="Select All",
                                method="restyle"
                            )
                        ]),
                        pad={"r": 10, "t": 10},
                        showactive=False,
                        x=1,
                        xanchor="right",
                        y=1.1,
                        yanchor="top"
                    ),
                ]
            )
        )

        fig.update_layout(
            xaxis_title=f"{level}",
            yaxis_title="Frequency"
        )

        # Show the plot
        figures.append(fig)



    # ------------------------------------------------------------------------------------------------------------------
    # Generate HTML ----------------------------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------------------------------

    with open('rna_qc_donor_mqc.html', 'w') as f:
        for fig in figures:
            f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))


print("[FINISH]")
# %%
