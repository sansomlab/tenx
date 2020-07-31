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
import anndata
import scanpy as sc
import pandas as pd
from scipy import sparse
import logging
import sys



# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

#logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger(__name__)
log_handler = logging.StreamHandler(sys.stdout)
log_handler.setFormatter(logging.Formatter('%(asctime)s %(message)s'))
log_handler.setLevel(logging.INFO)
L.addHandler(log_handler)
L.setLevel(logging.INFO)

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()


# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--anndata", default="anndata.h5ad", type=str,
                    help="File with the cell barcodes")
parser.add_argument("--mindist",default=0.5, type=str,
                    help="the minimum distance for the umap")
parser.add_argument("--outdir",default=1, type=str,
                    help="path to output directory")


# parser.add_argument("--resolution", default=1, type=str,
#                     help="the clustering resolution")

args = parser.parse_args()

L.info("Running with arguments:")
print(args)

# ########################################################################### #
# ############## Create outdir and set results file ######################### #
# ########################################################################### #


# write folder

# figures folder
sc.settings.figdir = args.outdir

sc.settings.set_figure_params(dpi=300, dpi_save=300)

# ########################################################################### #
# ############################### Run PAGA ################################## #
# ########################################################################### #

adata = anndata.read_h5ad(args.anndata)

# compute clusters
sc.tl.umap(adata,
           min_dist = float(args.mindist))

umap  = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP_1", "UMAP_2"])
umap["barcode"] = adata.obs["barcode"].values


umap.to_csv(os.path.join(args.outdir,
                      "umap." + args.mindist + ".tsv.gz"), sep="\t",
         index=False)


L.info("Complete")
