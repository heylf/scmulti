#!/usr/bin/env python

import argparse
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import random
import snapatac2 as snap
import anndata as ad
import numpy as np
import pandas as pd

# pip install nbformat
# pip install snapatac2
# pip install plotly
# pip install muon

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
parser.add_argument("-f", "--files", dest="files", metavar='str', required=True, help="Output dir")
parser.add_argument("-g", "--genome", dest="genome", metavar='str', required=True, help="Reference genome that was \
                    used for alignment.")

# Optional arguments
parser.add_argument("--demux", dest="demux", metavar='str', required=False, help="Demultiplexed files")

args = vars(parser.parse_args())

if __name__ == '__main__':

    print("[START]")

    # Seeds (!!! DO NOT CHANGE THIS SEED !!!)
    seed = 123
    random.seed(seed)
    print(f"[NOTE] seed {seed}")

    # Set genome
    genome = None

    if ( args['genome'] == "hg38" ):
        genome = snap.genome.hg38
    elif ( args['genome'] == "hg19" ):
        genome = snap.genome.hg19
    elif ( args['genome'] == "mm10" ):
        genome = snap.genome.mm10

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

    ####################################################################################################################
    ###### LOAD DATA ###################################################################################################
    ####################################################################################################################
    # This will get us the data from cellranger-arc
    # It also generates a table for the necessary information (arc_aggr_df.to_csv)
    print("[TASK] Load data")

    samples = args['samples'].split(",")

    fragment_files = []
   
    for file in args['files'].split(","):
        fragment_files.append("/".join(file.split("/")[:-1]) + "/atac_fragments.tsv.gz")

    atacs = snap.pp.import_data(fragment_files, 
                                chrom_sizes=genome,
                                sorted_by_barcode=False)

    # TODO turn histograms into lines (KDE)
    ####################################################################################################################
    ###### QC PLOTS ATAC ###############################################################################################
    ####################################################################################################################
    print("[TASK] Generate QC plots RNA")

    figures = []

    import plotly.subplots as sp

    combined_fig = sp.make_subplots(rows=1, cols=1)

    for i, atac in enumerate(atacs):
        fig = snap.pl.frag_size_distr(atac, show=False, interactive=True)
        fig.data[-1].name = samples[i]
        combined_fig.add_trace(fig.data[0], row=1, col=1)

    combined_fig.update_layout(title="Fragment Size Distributions")
    combined_fig.update_yaxes(title="Count")
    combined_fig.update_xaxes(title="Fragment size")
    figures.append(combined_fig)

    combined_fig = sp.make_subplots(rows=1, cols=1)

    for atac in atacs:
        fig = snap.pl.frag_size_distr(atac, show=False, interactive=True)
        fig.data[-1].name = samples[i]
        combined_fig.add_trace(fig.data[0], row=1, col=1)

    combined_fig.update_layout(title="log Fragment Size Distributions")
    combined_fig.update_yaxes(type="log", title="Count")
    combined_fig.update_xaxes(title="log Fragment size")
    figures.append(combined_fig)

    snap.metrics.tsse(atacs, genome)
    
    combined_fig = sp.make_subplots(rows=1, cols=1)  

    for i, atac in enumerate(atacs):
        fig = snap.pl.tsse(atac, interactive=True, show=False, out_file=None)
        fig.data[-1].name = samples[i]
        fig.data[0].showlegend = True
        combined_fig.add_trace(fig.data[0], row=1, col=1)

    # TODO bug in plot. You have to separate each of those plot and you cannot combine them!
    combined_fig.update_layout(title=f"TSS enrichment and number of unique fragments", 
                               legend=dict(orientation="h", y=1.02, x=0.5))
    combined_fig.update_yaxes(title="TSS enrichment score")
    combined_fig.update_xaxes(type="log", title="log Number of unique fragments")
    figures.append(combined_fig)

    for i, atac in enumerate(atacs):
        atac.obs["sample"] = samples[i]
        atac.obs.index = atac.obs.index + "_" + samples[i]
        atacs[i] = atac

    combined_df = pd.concat([atac.obs for atac in atacs])

    combined_df['log_n_fragment'] = np.log10(combined_df['n_fragment'] + 1)

    for level in ["n_fragment", "tsse"]:

        fig = px.histogram(combined_df, x=level, nbins=100,
                           labels={level: f"{level}"},
                           width=800, height=800)
        fig.update_layout(title=f"Distribution of {level} for all samples",
                          showlegend=True)
        figures.append(fig)

        if ( level != "tsse"):

            # Create a histogram trace for each sample
            traces = []
            for sample in samples:
                sample_data = combined_df[combined_df['sample'] == sample]
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
            sample_data = combined_df[combined_df['sample'] == sample]
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


# In[]
# ----------------------------------------------------------------------------------------------------------------------
# Generate HTML --------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

with open('2_atac_qc_sample_mqc.html', 'w') as f:
    for fig in figures:
        f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))


# In[]
########################################################################################################################
###### QC PLOTS ATAC DEMULIPLEXED ######################################################################################
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

    demux_combined_df = assignments.join(combined_df, how="left")


    # Get number of donors
    ndonor = len(set(demux_combined_df["donor_id"].tolist()))

    # ------------------------------------------------------------------------------------------------------------------
    # Plots ------------------------------------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------------------------------

    # Get unique donors
    unique_donors = demux_combined_df['donor_id'].unique()

    figures = []

    for level in ["n_fragment", "tsse"]:

        if ( level != "tsse"):

            # Create a histogram trace for each sample
            traces = []
            for donor in unique_donors:
                donor_data = demux_combined_df[demux_combined_df['donor_id'] == donor]
                trace = go.Histogram(x=donor_data[f'log_{level}'], nbinsx=100, name=donor)
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
                xaxis_title=f"log_{level}",
                yaxis_title="Frequency"
            )

            figures.append(fig)

        # Create a histogram trace for each sample
        traces = []
        for donor in unique_donors:
            donor_data = demux_combined_df[demux_combined_df['donor_id'] == donor]
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

    # ------------------------------------------------------------------------------------------------------------------
    # Generate HTML ----------------------------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------------------------------

    with open('4_atac_qc_donor_mqc.html', 'w') as f:
        for fig in figures:
            f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))


print("[FINISH]")
# %%
