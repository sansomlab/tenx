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

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger("run_paga")


sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()


# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #



parser = argparse.ArgumentParser()
parser.add_argument("--reduced_dims_matrix_file", default="reduced_dims.tsv.gz", type=str,
                    help="File with reduced dimensions")
parser.add_argument("--barcode_file", default="barcodes.tsv.gz", type=str,
                    help="File with the cell barcodes")
parser.add_argument("--outdir",default=1, type=str,
                    help="path to output directory")
parser.add_argument("--cluster_assignments", default=1, type=str,
                    help="gzipped tsv file with cell cluster assignments")
parser.add_argument("--comps", default="1", type=str,
                    help="Number of dimensions to include in knn and umap computation")
parser.add_argument("--k", default=20, type=int,
                    help="number of neighbors")

args = parser.parse_args()







# ########################################################################### #
# ############## Create outdir and set results file ######################### #
# ########################################################################### #


# write folder
results_file = args.outdir + "/" + "paga_anndata.h5ad"

# figures folder
sc.settings.figdir = args.outdir

sc.settings.set_figure_params(dpi=300, dpi_save=300)

# ########################################################################### #
# ############################### Run PAGA ################################## #
# ########################################################################### #


# Read matrix of reduced dimensions, create anndata and add dimensions
reduced_dims_mat = pd.read_csv(args.reduced_dims_matrix_file, sep="\t")
reduced_dims_mat.index = [x for x in pd.read_csv(args.barcode_file, header=None)[0]]

adata = sc.AnnData(obs=[x for x in reduced_dims_mat.index])
adata.obs.index = reduced_dims_mat.index
adata.obs.rename(columns={0:'barcode'}, inplace=True)

# Add dimensions to anndata
colnames = list(reduced_dims_mat.columns)
colname_prefix = list(set([re.sub(r'_[0-9]+', '', i) for i in colnames]))[0] + "_"
select_comps = args.comps.split(",")
select_comps = [ colname_prefix + item for item in select_comps]
reduced_dims_mat = reduced_dims_mat[select_comps]

L.info("Using comps " + ', '.join(list(reduced_dims_mat.columns)))

adata.obsm['X_pca'] = reduced_dims_mat.to_numpy(dtype="float32")

# Read and add cluster ids
df = pd.read_csv(args.cluster_assignments,sep="\t")
df.index = df["barcode"]

# Ensure correct ordering
adata.obs['cluster_id'] = df.loc[adata.obs.index,"cluster_id"].astype("category").values

# Run neighbors
L.info( "Using " + str(args.k) + " neighbors")

sc.pp.neighbors(adata, n_neighbors = args.k, use_rep = 'X_pca')



# Run, plot and save umap
sc.tl.umap(adata)

result  = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP_1", "UMAP_2"])
result["barcode"] = adata.obs["barcode"]

result.to_csv(os.path.join(args.outdir,
                           "umap.txt.gz"))

L.info("Complete")
