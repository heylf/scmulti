#!/usr/bin/env python

# In[]

import argparse
import muon as mu
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

# In[]
#########################################
###### GLOBAL VARS and DIRECTORIES ######
#########################################

# tool_description = """
# Runnin QC
# """

# # parse command line arguments
# parser = argparse.ArgumentParser(description=tool_description, formatter_class=argparse.RawDescriptionHelpFormatter)

# # version
# parser.add_argument("-o", "--out", dest="out", metavar='str', required=True, help="Single cell files")
# parser.add_argument("-v", "--version", action="version", version="%(prog)s 0.1.0")
# parser.add_argument("-s", "--samples", dest="samples", metavar='str', required=True, help="Sample names")
# parser.add_argument("-f", "--files", dest="files", metavar='str', required=True, help="Output dir")
# parser.add_argument("-g", "--genome", dest="genome", metavar='str', required=True, help="Reference genome that was \
# used for alignment. \")

# args = vars(parser.parse_args())

if __name__ == '__main__':

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


    # samples = args['samples'].split(",")
    samples = ['test_scARC_1', 'test_scARC_2']

    num_samples = 0

    fragment_files = []
    #for file in args['files'].split(","):
    for file in ["/home/florian/Documents/tmp_data_folder/output/cellrangerarc/count/test_scARC/outs/filtered_feature_bc_matrix.h5", "/home/florian/Documents/tmp_data_folder/output/cellrangerarc/count/test_scARC/outs/filtered_feature_bc_matrix.h5"]:
        fragment_files.append("/".join(file.split("/")[:-1]) + "/atac_fragments.tsv.gz")

    atacs = snap.pp.import_data(fragment_files, 
                                chrom_sizes=snap.genome.mm10,
                                sorted_by_barcode=False)

    # In[]
    # Test
    fragment_file = snap.datasets.pbmc5k()
    atacs = snap.pp.import_data([fragment_file, fragment_file], 
                                chrom_sizes=snap.genome.hg38,
                                sorted_by_barcode=False)
    samples = ['test_scARC_1', 'test_scARC_2']


    # In[]
    ###########################
    ###### QC PLOTS ATAC ######
    ###########################
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

    snap.metrics.tsse(atacs, snap.genome.mm10)
    
    combined_fig = sp.make_subplots(rows=1, cols=1)  

    for i, atac in enumerate(atacs):
        fig = snap.pl.tsse(atac, interactive=True, show=False, out_file=None)
        fig.data[-1].name = samples[i]
        fig.data[0].showlegend = True
        combined_fig.add_trace(fig.data[0], row=1, col=1)

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

        # TODO make more beautiful
        if ( level != "tsse"):
            fig = px.histogram(combined_df, x=f'log_{level}', nbins=100,
                            facet_col="sample", facet_col_wrap=4,
                            log_x=True, labels={level: f"log10 log_{level}"},
                            width=800, height=800)
            fig.update_layout(title=f"Distribution of log10 {level}",
                            showlegend=True)
            fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
            figures.append(fig)

        fig = px.histogram(combined_df, x=level, nbins=100,
                           facet_col="sample", facet_col_wrap=4,
                           log_x=False, labels={level: f"{level}"},
                           width=800, height=800)
        fig.update_layout(title=f"Distribution of {level}",
                          showlegend=True)
        figures.append(fig)
        fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
        figures.append(fig)

    # In[]
    # ----------------------------------------------------------------------------------------------------------------------
    # Generate HTML --------------------------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------------

    with open('atac_qc_mqc.html', 'w') as f:
        for fig in figures:
            f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))

    print("[FINISH]")
    # %%
