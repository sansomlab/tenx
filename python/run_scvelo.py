import os
import logging
import sys
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')

from scipy import io
import scanpy as sc
import anndata as ad

import scvelo as scv


import argparse

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
sc.logging.print_versions()


# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--loom", default="none", type=str,
                    help="A loom file")
parser.add_argument("--dropest_dir", default="none", type=str,
                    help="File with dropest layers")
parser.add_argument("--outdir",default=1, type=str,
                    help="path to output directory")
parser.add_argument("--cluster_assignments", default=1, type=str,
                    help="gzipped tsv file with cell cluster assignments")
parser.add_argument("--cluster_colors", default=1, type=str,
                    help="tsv file with the color palette for the clusters")
parser.add_argument("--resolution", default=1, type=str,
                    help="the clustering resolution")
parser.add_argument("--rdims", default="1", type=str,
                    help="reduced dimensions file")
parser.add_argument("--rdim_method", default="umap", type=str,
                    help="name of the dimension reduction method")
parser.add_argument("--rdim1", default="UMAP1", type=str,
                    help="reduced dimension 1")
parser.add_argument("--rdim2", default="UMAP2", type=str,
                    help="reduced dimension 2")
parser.add_argument("--barcodes_dir", default=1, type=str,
                    help="gzipped tsv file with barcodes to keep")

args = parser.parse_args()

# ########################################################################### #
# ######################## Initialise AnnData ############################### #
# ########################################################################### #

if not args.loom == "none":

    adata = scv.read(args.loom)
    # get directory with metadata + barcodes
    metadata_dir = args.rdims.split("/")[0]

elif not args.dropest_dir == "none":

    exon_matrix = os.path.join(args.dropest_dir, "exons.mtx.gz")
    intron_matrix = os.path.join(args.dropest_dir, "introns.mtx.gz")
    spanning_matrix =  os.path.join(args.dropest_dir, "spanning.mtx.gz")

    exons = io.mmread(exon_matrix).transpose().tocsr()
    introns = io.mmread(intron_matrix).transpose().tocsr()
    spanning = io.mmread(spanning_matrix).transpose().tocsr()

    adata = ad.AnnData(X=exons)
    adata.layers["spliced"] = adata.X
    adata.layers["unspliced"] = introns
    adata.layers["ambiguous"] = spanning

    adata.obs.index = [x for x in
                       pd.read_csv(os.path.join(args.dropest_dir, "barcodes.tsv.gz"),
                                   header=None)[0].values]
    metadata_dir = args.dropest_dir

else:
    raise ValueError("either a loom file or dropEst directory must be specified")

# Add the variable and observation information

feat = pd.read_csv(os.path.join(metadata_dir,"features.tsv.gz"), header=None)
feat.columns = ["gene_symbol"]
feat.index = feat["gene_symbol"]
samples  = pd.read_csv(os.path.join(metadata_dir, "barcodes.tsv.gz"),header=None)
metadata = pd.read_csv(os.path.join(metadata_dir,"metadata.tsv.gz"), sep="\t")
metadata.index = metadata.barcode.values

adata.vars = feat
adata.obs = metadata.loc[adata.obs.index,]

# read barcodes_to_keep (from main subsetted matrix)
barcodes_to_keep = pd.read_csv(os.path.join(args.barcodes_dir, "barcodes.tsv.gz"),header=None)

# convert dataframe to array
tmp = barcodes_to_keep.to_numpy()
tmpf = tmp.flatten()

# subset adata
adata = adata[tmpf]

L.info("AnnData initialised")

# ########################################################################### #
# ################### Get rdims (e.g. umap) info  ########################### #
# ########################################################################### #


# ########################################################################### #
# ############################ run scvelo ################################### #
# ########################################################################### #

scv.settings.figdir = args.outdir
print(scv.settings.figdir)

# Get the color palette
# ggplot_palette = [x for x in pd.read_csv(args.cluster_colors,
#                       header=None, sep="\t")[0].values]

# preprocessing
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)

# velocity computation
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)

L.info("Velocity calculation completed")


cluster_files = [x.strip() for x in args.cluster_assignments.split(",")]
resolutions = [x.strip() for x in args.resolution.split(",")]
color_files =  [x.strip() for x in args.cluster_colors.split(",")]
rdims_files =  [x.strip() for x in args.rdims.split(",")]

clusterings = dict(zip(resolutions, cluster_files))
colorings = dict(zip(resolutions, color_files))
rdims_tables = dict(zip(resolutions, rdims_files * len(resolutions)))

for resolution, cluster_file in clusterings.items():

    print(rdims_tables.keys())
    print(colorings.keys())
    rdims = pd.read_csv(rdims_tables[resolution],
                        sep="\t")

    rdims.index = rdims["barcode"]

    # adata.obs["clusters"] = rdims["cluster"].values
    adata.obsm["X_" + args.rdim_method] = rdims.loc[adata.obs.index,
                                                [args.rdim1, args.rdim2]].values

    L.info("Rdims info added")


    clusters = pd.read_csv(cluster_file,sep="\t")

    # Get the color palette
    ggplot_palette = [x for x in pd.read_csv(colorings[resolution],
                      header=None, sep="\t")[0].values]

    # ggplot_cmap = ListedColormap(sns.color_palette(ggplot_palette).as_hex())

    clusters = pd.read_csv(cluster_file, sep="\t")
    clusters.index = [x for x in clusters.barcode.values]

    adata = adata[clusters.index]

    adata.obs['cluster'] = clusters.loc[adata.obs.index,
                                        "cluster_id"].astype("category").values

    # make the plots
    fname = args.rdim_method + "_stream." + resolution + ".png"
    scv.pl.velocity_embedding_stream(adata,
                                     basis=args.rdim_method,
                                     dpi=300,
                                     color="cluster",
                                     palette=ggplot_palette,
                                     save=fname,
                                     show=False)

    fname = args.rdim_method + "_velocity." + str(resolution) + ".png"
    scv.pl.velocity_embedding(adata, basis=args.rdim_method,
                              arrow_length=2, arrow_size=1.5,
                              dpi=300,
                              color="cluster", palette=ggplot_palette,
                              save=fname,
                              show=False)

    fname = args.rdim_method + "_grid." + str(resolution) + ".png"
    scv.pl.velocity_embedding_grid(adata, basis=args.rdim_method,
                                   dpi=300,
                                   color="cluster", palette=ggplot_palette,
                                   save=fname,
                                   show=False)


L.info("Plotting finished")
