import os
import re
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')

from matplotlib import rcParams
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as pl
import seaborn as sns
import scanpy as sc
import pandas as pd
from scipy import sparse
import logging
import sys



# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

L = logging.getLogger(__name__)
log_handler = logging.StreamHandler(sys.stdout)
log_handler.setFormatter(logging.Formatter('%(asctime)s %(message)s'))
log_handler.setLevel(logging.INFO)
L.addHandler(log_handler)
L.setLevel(logging.INFO)

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
# sc.logging.print_versions()


# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #



parser = argparse.ArgumentParser()
# parser.add_argument("--reduced_dims_matrix_file", default="reduced_dims.tsv.gz", type=str,
#                     help="File with reduced dimensions")
# parser.add_argument("--barcode_file", default="barcodes.tsv.gz", type=str,
#                     help="File with the cell barcodes")
parser.add_argument("--anndata", default="anndata.h5ad",
                    help="The anndata object")
# parser.add_argument("--umap", default="umap.tsv.gz",
#                    help="The umap coordinates")
# parser.add_argument("--features_file", default="features.tsv.gz", type=str,
#                     help="File with the feature names")
parser.add_argument("--outdir",default=1, type=str,
                    help="path to output directory")
parser.add_argument("--cluster_assignments", default=1, type=str,
                    help="gzipped tsv file with cell cluster assignments")
parser.add_argument("--cluster_colors", default=1, type=str,
                    help="tsv file with the color palette for the clusters")
# parser.add_argument("--comps", default="1", type=str,
#                     help="Number of dimensions to include in knn and umap computation")
# parser.add_argument("--k", default=20, type=int,
#                     help="number of neighbors")

args = parser.parse_args()







# ########################################################################### #
# ############## Create outdir and set results file ######################### #
# ########################################################################### #


# write folder
# results_file = args.outdir + "/" + "paga_anndata.h5ad"

# figures folder
sc.settings.figdir = args.outdir

# Get the color palette
ggplot_palette = [x for x in pd.read_csv(args.cluster_colors,
                      header=None, sep="\t")[0].values]

ggplot_cmap = ListedColormap(sns.color_palette(ggplot_palette).as_hex())

sc.settings.set_figure_params(dpi=300, dpi_save=300)

# ########################################################################### #
# ############################### Run PAGA ################################## #
# ########################################################################### #

# Read in the anndata object with pre-computed neighbors.
adata = sc.read_h5ad(args.anndata)

# Read and add cluster ids
df = pd.read_csv(args.cluster_assignments,sep="\t")
df.index = df["barcode"]

# Ensure correct ordering
adata.obs['cluster_id'] = df.loc[adata.obs.index,"cluster_id"].astype("category").values

print(ggplot_palette)

# Run and plot paga
sc.tl.paga(adata, groups='cluster_id')
sc.pl.paga(adata, save=".png", show=False, cmap=ggplot_cmap)

# Run, plot and store paga-initialised umap
sc.tl.umap(adata, init_pos = 'paga')

sc.pl.umap(adata, color="cluster_id", legend_loc='on data',
           save = ".paga.initialised.png", show=False,
           palette=ggplot_palette)

# Save paga-initialised UMAP coordinates
umap  = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP_1", "UMAP_2"])
umap["barcode"] = adata.obs["barcode"].values

umap.to_csv(os.path.join(args.outdir,
                         "umap.paga.init.tsv.gz"), sep="\t",
            index=False)

# Compute and plot the force directed graph (FDG)
sc.tl.draw_graph(adata)
sc.pl.draw_graph(adata, color='cluster_id', legend_loc='on data',
                 save=".png", show=False, palette=ggplot_palette)

# Compute and plot the PAGA initialised FDG
sc.tl.draw_graph(adata, init_pos='paga')
sc.pl.draw_graph(adata, color='cluster_id', legend_loc='on data',
                 save=".paga.initialised.png", show=False,
                 palette=ggplot_palette)

paga_fa2 = pd.DataFrame(adata.obsm["X_draw_graph_fa"],
                                    columns=["FA1","FA2"])

paga_fa2["barcode"] = adata.obs["barcode"].values

paga_fa2.to_csv(os.path.join(args.outdir,
                             "paga_init_fa2.tsv.gz"),
                             sep="\t")

L.info("Complete")
