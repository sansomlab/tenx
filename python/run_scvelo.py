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
parser.add_argument("--dropest_dir", default=1, type=str,
                    help="File with dropest layers")
parser.add_argument("--outdir",default=1, type=str,
                    help="path to output directory")
parser.add_argument("--cluster_assignments", default=1, type=str,
                                        help="gzipped tsv file with cell cluster assignments")
parser.add_argument("--cluster_colors", default=1, type=str,
                    help="tsv file with the color palette for the clusters")
parser.add_argument("--rdims", default="1", type=str,
                    help="reduced dimensions file")
parser.add_argument("--rdim_method", default="umap", type=str,
                    help="name of the dimension reduction method")
parser.add_argument("--rdim1", default="UMAP1", type=str,
                    help="reduced dimension 1")
parser.add_argument("--rdim2", default="UMAP2", type=str,
                    help="reduced dimension 2")


args = parser.parse_args()

# ########################################################################### #
# ######################## Initialise AnnData ############################### #
# ########################################################################### #

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

# Add the variable and observation information
feat = pd.read_csv(os.path.join(args.dropest_dir,"features.tsv.gz"), header=None)
feat.columns = ["gene_symbol"]
feat.index = feat["gene_symbol"]
samples  = pd.read_csv(os.path.join(args.dropest_dir, "barcodes.tsv.gz"),header=None)
metadata = pd.read_csv(os.path.join(args.dropest_dir,"metadata.tsv.gz"), sep="\t")
metadata.index = metadata.barcode.values

adata.vars = feat
adata.obs = metadata

L.info("AnnData initialised")

clusters = pd.read_csv(args.cluster_assignments, sep="\t")
clusters.index = clusters["barcodes"]

adata.obs['cluster'] = clusters.loc[adata.obs.index,"cluster_id"].astype("category").values

# ########################################################################### #
# ################### Get rdims (e.g. umap) info  ########################### #
# ########################################################################### #

rdims = pd.read_csv(args.rdims, sep="\t")

rdims.index = rdims["barcode"]

rdata = adata[rdims.index]
# rdata.obs["clusters"] = rdims["cluster"].values
rdata.obsm["X_" + args.rdim_method] = rdims[[args.rdim1, args.rdim2]].values

L.info("Rdims info added")

# ########################################################################### #
# ############################ run scvelo ################################### #
# ########################################################################### #

scv.settings.figdir = args.outdir + "/"

# Get the color palette
ggplot_palette = [x for x in pd.read_csv(args.cluster_colors,
                      header=None, sep="\t")[0].values]

# preprocessing
scv.pp.filter_and_normalize(rdata)

scv.pp.moments(rdata)

# velocity computation
scv.tl.velocity(rdata)

scv.tl.velocity_graph(rdata)

L.info("Velocity calculation completed")

# make the plots
scv.pl.velocity_embedding_stream(rdata, basis=args.rdim_method,dpi=300,
                                 color="cluster", palette=ggplot_palette,
                                 save=args.rdim_method + "_stream.png",
                                 show=False)

scv.pl.velocity_embedding(rdata, basis=args.rdim_method,
                          arrow_length=2, arrow_size=1.5,
                          dpi=300,
                          color="cluster", palette=ggplot_palette,
                          save=args.rdim_method + "_velocity.png",
                          show=False)

scv.pl.velocity_embedding_grid(rdata, basis=args.rdim_method,
                               dpi=300,
                               color="cluster", palette=ggplot_palette,
                               save=args.rdim_method + "_grid.png",
                               show=False)


L.info("Plotting finished")
