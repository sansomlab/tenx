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
import pandas as pd
from scipy import sparse
import logging
import sys
import scprep
import phate
import anndata


# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

L = logging.getLogger(__name__)
log_handler = logging.StreamHandler(sys.stdout)
log_handler.setFormatter(logging.Formatter('%(asctime)s %(message)s'))
log_handler.setLevel(logging.INFO)
L.addHandler(log_handler)
L.setLevel(logging.INFO)


# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--data", default="data.tsv.gz", type=str,
                    help="File with the data, e.g. scale.data from Seurat")
parser.add_argument("--assay", default="reduced.dimensions", type=str,
                    help="type of assay: scaled.data or dimension.reduction")
parser.add_argument("--barcode_file", default="barcodes.tsv.gz", type=str,
                    help="File with the cell barcodes")
parser.add_argument("--outdir",default=1, type=str,
                    help="path to output directory")
parser.add_argument("--resolution",default=1, type=str,
                    help="comma separated list of cluster resolutions")
parser.add_argument("--cluster_assignments", default=1, type=str,
                    help="comma separated list of gzipped tsv files with cell cluster assignments")
parser.add_argument("--cluster_colors", default=1, type=str,
                    help="tsv file with the color palette for the clusters")
parser.add_argument("--k", default=5, type=int,
                    help="number of neighbors")
parser.add_argument("--gif", default="No", type=str,
                    help="output a GIF")
parser.add_argument("--reduction_name", default="pca", type=str,
                    help="Name of the reduced dimseions")
parser.add_argument("--input_type", default="tsv", type=str,
                    help="anndata or tsv")

args = parser.parse_args()

L.info("Running with arguments:")
print(args)


# ########################################################################### #
# ############################## Run PHATE ################################## #
# ########################################################################### #

if args.input_type == "tsv":
    if args.assay == "reduced.dimensions":
        # Read matrix of reduced dimensions, create anndata and add dimensions
        data = pd.read_csv(args.data, sep="\t", header=0)

    if args.assay == "scaled.data":
        # Read matrix of reduced dimensions, create anndata and add dimensions
        data = pd.read_csv(args.data, sep="\t", header=None)
        # we need to transpose the data for PHATE
        data = data.transpose()

if args.input_type == "anndata":
    # Read in anndata and extract correct data type
    adata = anndata.read(args.data)
    if args.assay == "reduced.dimensions":
        obsm_use = 'X_' + str(args.reduction_name)
        data = pd.DataFrame(adata.obsm[obsm_use].copy(),
                            index=adata.obs.index.copy())
    if args.assay == "scaled.data":
        data = pd.DataFrame(adata.X.copy(),
                            index=adata.index.copy())
        data = data.transpose()

phate_operator = phate.PHATE(n_jobs=-2, knn=args.k)
x2 = phate_operator.fit_transform(data)

# Read and add cluster ids
cluster_files = [x.strip() for x in args.cluster_assignments.split(",")]
resolutions = [x.strip() for x in args.resolution.split(",")]
color_files =  [x.strip() for x in args.cluster_colors.split(",")]

clusterings = dict(zip(resolutions, cluster_files))
colorings = dict(zip(resolutions, color_files))

for resolution, cluster_file in clusterings.items():

    clusters = pd.read_csv(cluster_file,sep="\t")

    # Get the color palette
    ggplot_palette = [x for x in pd.read_csv(colorings[resolution],
                      header=None, sep="\t")[0].values]

    ggplot_cmap = ListedColormap(sns.color_palette(ggplot_palette).as_hex())


    # save a 2D plot
    scprep.plot.scatter2d(x2, c=clusters["cluster_id"],
                      figsize=(12,8), cmap=ggplot_cmap,
                      ticks=False, label_prefix="PHATE", s=15,
                      filename=os.path.join(args.outdir,"phate." + resolution + ".2D.png"),
                      dpi=300)

# save the 2D coordinates
rdims_phate = pd.DataFrame(x2,
                           columns=["PHATE1","PHATE2"])

if args.input_type == "tsv":
    rdims_phate["barcode"] = pd.read_csv(args.barcode_file, header=None)[0].values
if args.input_type == "anndata":
    rdims_phate["barcode"] = adata.obs.index.values


rdims_phate.to_csv(os.path.join(args.outdir,"phate.tsv.gz"),
                   sep="\t")

# save a 3D plot
phate_operator.set_params(n_components=3)
x3 = phate_operator.transform()


for resolution, cluster_file in clusterings.items():

    clusters = pd.read_csv(cluster_file,sep="\t")

    # Get the color palette
    ggplot_palette = [x for x in pd.read_csv(colorings[resolution],
                      header=None, sep="\t")[0].values]

    ggplot_cmap = ListedColormap(sns.color_palette(ggplot_palette).as_hex())


    scprep.plot.scatter3d(x3, c=clusters["cluster_id"],
                          figsize=(8,6), cmap=ggplot_cmap,
                          ticks=False, label_prefix="PHATE",
                          filename=os.path.join(args.outdir,"phate." + resolution + ".3D.png"),
                          dpi=300)


# save the 3D coordinates
rdims_phate = pd.DataFrame(x3,
                           columns=["PHATE1","PHATE2", "PHATE3"])

if args.input_type == "tsv":
    rdims_phate["barcode"] = pd.read_csv(args.barcode_file, header=None)[0].values
if args.input_type == "anndata":
    rdims_phate["barcode"] = adata.obs.index.values

rdims_phate.to_csv(os.path.join(args.outdir,"phate_3D.tsv.gz"),
                   sep="\t")

# save a GIF!
if args.gif.lower() == "yes":

    for resolution, cluster_file in clusterings.items():

        clusters = pd.read_csv(cluster_file,sep="\t")

        # Get the color palette
        ggplot_palette = [x for x in pd.read_csv(colorings[resolution],
                      header=None, sep="\t")[0].values]

        ggplot_cmap = ListedColormap(sns.color_palette(ggplot_palette).as_hex())


        scprep.plot.rotate_scatter3d(x3, c=clusters["cluster_id"],
                                 figsize=(8,6), cmap=ggplot_cmap,
                                 ticks=False, label_prefix="PHATE",
                                 filename=os.path.join(args.outdir,"phate." + resolution + ".3D.gif"),
                                 dpi=300)


L.info("Phate accompli")
